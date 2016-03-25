/*
 * kmpare.cpp
 */

#include <iostream>
#include <cstring>
#include <cstdlib>
#include <sstream>
#include "kmpare.h"
#include "kmer.h"
#include "parseData.h"
#include "memPool.h"

//int main (int argc, char** argv)
int main (int argc, char** argv)
{
	// initialize variables and get user input
	std::vector<std::string> infiles;
	std::vector< std::vector<unsigned int> > sets;
	std::string fout;
	if (argc == 1)
	{
		info(version);
		return 0;
	}
	if ( !parseArgs(argc, argv, &infiles, &sets, fout) )
		return 0;

	// initialize objects
	kmer jellydata; // handles kmer data

	// open outfile stream
	if ( fexists(fout.c_str()) )
	{
		std::cerr << "File already exists: " << fout << "\n" << "-->exiting";
		return 0;
	}
	std::ofstream os(fout.c_str());
	if (os.fail())
	{
		std::cerr << "Could not open file: " << fout << "\n" << "-->exiting\n";
		return 1;
	}

	// parse Jellyfish files
	jellydata.parseJellyCounts(infiles);
	if (jellydata.fail)
	{
		std::cerr << "--> exiting\n";
		return 1;
	}
	std::cerr << jellydata.nkmers() << " kmer sequences in the dataset\n";

	// analyze kmer counts
	MemPool<double> stats;
	jellydata.fit(&jellydata.datamap, &stats, &sets);
	if (jellydata.fail)
	{
		std::cerr << "ERROR: Kmer count analysis failed\n" << "--> exiting\n";
		return 1;
	}

	// print result
	std::cerr << "Dumping results to file: " << fout << "\n";
	printHeader(os, infiles.size(), &sets);
	jellydata.printStats(os, &jellydata.datamap, sets.size(), &stats);
	if (jellydata.fail)
	{
		std::cerr << "ERROR: Printing results failed\n" << "--> exiting\n";
		return 1;
	}

	std::cerr << "finished!\n";
	return 0;
}

bool parseArgs (int argc, char** argv, std::vector<std::string>* ifname, std::vector< std::vector<unsigned int> >* cmpindex, std::string& ofname)
{
	int argpos = 1;
	int counter = 0;

	// get arguments
	while (argpos < argc)
	{
		if ( strcmp(argv[argpos], "-help") == 0 )
		{
			info(version);
			return false;
		}
		else if ( strcmp(argv[argpos], "-infile") == 0 )
		{
			++argpos;
			ifname->reserve(2);
			while ( strcmp(argv[argpos], "-compset") && strcmp(argv[argpos], "-outfile"))
			{
				ifname->push_back(argv[argpos]);
				++counter;
				++argpos;
				if (argpos >= argc)
					break;
			}
			ifname->resize(counter);
			counter = 0;
		}
		else if ( strcmp(argv[argpos], "-compset") == 0)
		{
			++argpos;
			cmpindex->reserve(1);
			while ( strcmp(argv[argpos], "-infile") && strcmp(argv[argpos], "-outfile"))
			{
				if (argv[argpos][0] == '{')
				{
					cmpindex->push_back(parseSet(argc, argv, argpos));
					++counter;
				}
				++argpos;
				if (argpos >= argc)
					break;
			}
			cmpindex->resize(counter);
			counter = 0;
		}
		else if ( strcmp(argv[argpos], "-outfile") == 0)
		{
			ofname = argv[argpos + 1];
			argpos += 2;
		}
		else
		{
			fprintf(stderr, "Unknown command: %s\n", argv[argpos]);
			return false;
		}
	}

	//check arguments
	if (ofname.empty())
	{
		fprintf(stderr, "Must supply -outfile\n");
		return false;
	}

	if (ifname->empty())
	{
		fprintf(stderr, "Must supply -infile\n");
		return false;
	}

	if (!cmpindex->empty())
	{
		unsigned int nfiles_index = ifname->size() - 1;
		std::vector<unsigned int>::const_iterator elemIter;
		for(std::vector< std::vector<unsigned int> >::const_iterator setIter = cmpindex->begin(); setIter != cmpindex->end(); ++setIter)
		{
			for(elemIter = (*setIter).begin(); elemIter != (*setIter).end(); ++elemIter)
			{
				if ( *elemIter > nfiles_index)
				{
					fprintf(stderr, "One of the libraries supplied to -compset exceeds number of input files\n");
					return false;
				}
			}
		}
	}

	return true;
}

std::vector<unsigned int> parseSet (int argc, char** argv, int& pos)
{
	std::vector<unsigned int> set;
	set.reserve(2);
	int n = 0;
	std::string str = argv[pos];
	if (str.length() > 1)
		str = str.substr(1, str.length() - 1);
	while (str[str.length()-1] != '}' && pos < argc)
	{
		if (str != "{")
		{
			set.push_back( static_cast<unsigned int>(atoi(str.c_str()) - 1) );
			++n;
		}
		++pos;
		str = argv[pos];
	}
	if (str.length() > 1)
	{
		str = str.substr(0, str.length() - 1);
		set.push_back( static_cast<unsigned int>(atoi(str.c_str()) - 1) );
		++n;
	}
	set.resize(n);
	return set;
}

void printHeader (std::ofstream& os, unsigned int nlibs, const std::vector< std::vector<unsigned int> >* sets)
{
	os << "kmer";
	for(unsigned int i = 1; i <= nlibs; ++i)
	{
		os << "\t" << "lib" << i;
	}
	if (!sets->empty())
	{
		std::vector<unsigned int>::const_iterator libIter;
		for(std::vector< std::vector<unsigned int> >::const_iterator setIter = sets->begin(); setIter != sets->end(); ++setIter)
		{
			os << "\t" << "{ ";
			for(libIter = (*setIter).begin(); libIter != (*setIter).end(); ++libIter)
			{
				os << *libIter + 1 << " ";
			}
			os << "}";
		}
	}
	os << "\n";
}

void info (const char* v)
{
	fprintf(stderr, "\nkmpare version %s\n", v);
	std::cerr << "\nInput:\n"
	<< "-infile FILE Jellyfish text files of kmer counts\n"
	<< "-compset {INT} set(s) of libraries to compare\n"
	<< "-outfile FILE output file name\n"
	<< "\nOutput:\n"
	<< "<kmer> <library count> <goodness-of-fit for library set>\n"
	<< "\n";
}
