/*
 * kmer.cpp
 */

#include <cstdio>
#include <limits>
#include <math.h>
#include <algorithm>
#include "kmer.h"
#include "parseData.h"
#include <iostream>

kmer::kmer ()
	: fail(0),
	  stat(0),
	  statsize(0),
	  nonseq_char(3),
	  xtra_reserve(0.50),
	  nlibs(0),
	  kmertypes(0),
	  storage(0)
{

}

kmer::~kmer ()
{
	if (stat != 0)
		clearStat();
}

// parseJellyCounts extracts kmers and counts from Jellyfish files (sets member "counts")
void kmer::parseJellyCounts (std::vector<std::string>& files)
{
	if ( files.empty() )
	{
		fprintf(stderr, "No files to parse in call to kmer::parseJellyCounts\n");
		fail = 1;
		return;
	}
	unsigned int lib = 0;
	unsigned long int filelen = 0;
	Value<double> seqdat;
	seqdat.count.setSize(files.size());
	Key<long int> seqID;
	int maxdigit = ndigit(std::numeric_limits<long int>::max());
	std::ifstream is;
	for (std::vector<std::string>::iterator fIter = files.begin(); fIter != files.end(); ++fIter)
	{
		std::cerr << "reading file: " << *fIter << "\n";
		is.open(fIter->c_str());
		if (is.fail())
		{
				std::cerr << "Could not open file: " << *fIter << "\n";
				fail = 1;
				return;
		}
		if (is.rdbuf()->in_avail() == 0)
		{
			std::cerr << "0 sequences found in file: " << *fIter << "\n";
			fail = 1;
			return;
		}
		int merlen = jellyMerLength(is);
		fprintf(stderr, "kmer length is %d\n", merlen);
		int seqparts = ceil(merlen / static_cast<float>(maxdigit));
		if (lib == 0)
		{
			seqID.id.setSize(seqparts);
			libtotal.setSize(files.size());
			filelen = estLines (is, merlen, nonseq_char);
			storage = filelen + filelen * xtra_reserve;
			seqID.tableSize = &storage;
			datamap.reserve(storage);
		}
		std::string line;
		std::vector<std::string> tokens;
		std::pair<countmap::iterator, bool> result;
		size_t nbuckets = 0;
		unsigned int i = 0;
		while (getline(is, line))
		{
			if (datamap.size() > datamap.bucket_count())
			{
				std::cerr << "Increasing reserve to accommodate additional kmers...\n";
				nbuckets = datamap.size() + datamap.size() * xtra_reserve;
				datamap.reserve(nbuckets);
				storage = nbuckets;
			}
			tokens = split(line, ' ');
			seqdat.count[lib] = (atoi(tokens[1].c_str()));
			libtotal[lib] += seqdat.count[lib];
			seqtonum(tokens[0], seqID.id, maxdigit);
			result = datamap.insert(countmap::value_type(seqID, seqdat));
			if (result.second == 1)
			{
				++kmertypes;
				for (i = 0; i < lib; ++i)
					result.first->second.count[i] = 0;
			}
                        result.first->second.count[lib] = seqdat.count[lib];
		}
		is.close();
		++lib;
	}

	/* start debug code
	std::stringstream debugss;
	for(countmap::iterator it = datamap.begin(); it != datamap.end(); ++it)
	{
		std::cerr << numtoseq(it->first.id, debugss) << "\t";
		for (unsigned int j = 0; j < it->second.count.size(); ++j)
		{
			std::cerr << "\t" << it->second.count[j];
		}
		std::cerr << "\n";
	}
	 end debug code */
}
// numtoseq converts a string of numbers to nucleotide letters
std::string kmer::numtoseq (const Array<long int>& num, std::stringstream& ss) const
{
	if (!ss.str().empty())
		ss.str(std::string());

	unsigned int i = 0;
	for (i = 0; i < num.size(); ++i)
		ss << num[i];
	std::string s( ss.str() );
	for(i = 0; i < s.length(); ++i)
	{
		if (s[i] == '1')
			s[i] = 'A';
		else if ( s[i] == '2')
			s[i] = 'C';
		else if ( s[i] == '3')
			s[i] = 'G';
		else if ( s[i] == '4')
			s[i] = 'T';
		else if ( s[i] == '5')
			s[i] = 'N';
		else
		{
			fprintf(stderr, "WARNING: Unknown value in sequence\n");
			s[i] = '?';
		}
	}

	return s;
}



// seqtonum converts a string to a set of integers
void kmer::seqtonum (const std::string& s, Array<long int>& num, const int maxdigit)
{
	std::string::const_iterator iter;
	int i = 0;
	int j = 0;
	char charstr [maxdigit];
	char nuc;
	for (iter = s.begin(); iter != s.end(); ++iter)
	{
		nuc = toupper(*iter);
		if (i < maxdigit)
		{
			if (nuc == 'A')
				charstr[i] = '1';
			else if (nuc == 'C')
				charstr[i] = '2';
			else if (nuc == 'G')
				charstr[i] = '3';
			else if (nuc == 'T')
				charstr[i] = '4';
			else if (nuc == 'N')
				charstr[i] = '5';
			else
			{
				fprintf(stderr, "WARNING: Unknown base found in sequence\n");
				charstr[i] = '6';
			}
			i++;
		}
		else
		{
			num[j] = atol(charstr);
			i = 0;
			++j;
			--iter;
		}
	}
	for (int k = i; k < maxdigit; ++k)
		charstr[k] = 0;
	if (i < maxdigit)
		num[j] = atol(charstr);
}

// ndigit determines the number of digits in a number
template <class T> int kmer::ndigit (T number)
{
	int len = 0;
	while (number)
    {
		++len;
		number /= 10;
     }
    return len;
}

// jellyMerLength determines length of kmers in Jellyfish file
int kmer::jellyMerLength (std::ifstream& is)
{
	if (is)
	{
		int currentpos = is.tellg();
		std::string line;
		while (line.empty())
		{
			getline(is, line);
		}
		is.seekg(currentpos);
		int kmerlength = (split(line, ' '))[0].length();
		return kmerlength;
	}
	else
	{
		fprintf(stderr, "ERROR: No instream to function jellyMerLength\n");
		return -1;
	}
}

//estLines approximates the number of lines in Jellyfish file
unsigned long int kmer::estLines (std::ifstream& is, int merlength, const int nonseq_n)
{
	if (merlength <= 0)
	{
		fprintf(stderr, "ERROR: Invalid kmer length in estLines function");
		return -1;
	}
	if (is)
	{
		unsigned long int currentpos = is.tellg();
		is.seekg(0, is.end);
		unsigned long int length = is.tellg();
		is.seekg(currentpos);
		return ceil((length/static_cast<double>((merlength + nonseq_n*sizeof(char)))));
	}
	else
	{
		fprintf(stderr, "ERROR: No instream to function estLines\n");
		return -1;
	}
}


// calculates the sum for a vector given a set of vector indices
template <class T> T kmer::arraySum (const Array<T>& v, std::vector<unsigned int>* index)
{
	T sum = 0;
	if (!index)
		for (unsigned int i = 0; i < v.size(); ++i)
			sum += v[i];
	else
	{
		for(std::vector<unsigned int>::const_iterator indit = index->begin(); indit != index->end(); ++indit)
		{
			if (*indit > v.size() - 1)
			{
				fprintf(stderr, "ERROR: Number of indices exceeds array size\n");
				fail = 1;
				return -1.0/0.0;
			}
			sum += v[*indit];
		}
	}
	return sum;
}

// gets probability that a randomly chosen kmer from a pool of kmers comes from a particular library
void kmer::libProbs (double p [], std::vector<unsigned int>* idx, Array<size_t>& lib_count)
{
	size_t total = arraySum(lib_count, idx);
	unsigned int i = 0;
	static std::vector<unsigned int>::const_iterator iter;
	for (iter = idx->begin(); iter != idx->end(); ++iter)
	{
		p[i] = static_cast<double>(lib_count[*iter])/total;
		++i;
	}
}

// calculates goodness-of-fit statistic using a weighted average as the expected value
template <class T> double kmer::calcWGOF (double p [], const Array<T>& obs, std::vector<unsigned int>* idx)
{
	double stat = 0.0;
	double expt = 0;
	T obs_total = arraySum(obs, idx);
	unsigned int i = 0;
	for (std::vector<unsigned int>::const_iterator iter = idx->begin(); iter != idx->end(); ++iter)
	{
		expt = p[i] * obs_total;
		if (expt != 0)
			stat += pow((obs[*iter] - expt), 2) / expt;
		else
		{
			fprintf(stderr, "WARNING: Division by zero in calcGOF\n");
			return 1.0/0.0;
		}
		++i;
	}
	return stat;
}

// calculates goodness-of-fit for a set of kmer counts and stores them in a MemPool<double> object
void kmer::fit(countmap* data, MemPool<double>* stats, std::vector< std::vector<unsigned int> >* set)
{
	unsigned int j = 0;

	double* p [set->size()]; // P(kmer comes from library i)
	for (j = 0; j < set->size(); ++j)
	{
		p[j] = new double[(*set)[j].size()];
		libProbs(p[j], &(*set)[j], libtotal);
	}

	// allocate storage for GOF statistics
	std::cerr << "Allocating space for goodness-of-fit statistics...\n";
	if (stats->initial() == 1)
	{
		std::cerr << "WARNING: memPool object supplied to kmer::fit already initialized --> clearing content\n";
		stats->deleteReserve();
	}
	stats->formatReserve((kmertypes * set->size()), 0);

	// assign GOF values to memPool buffer
	std::cerr << "Calculating goodness-of-fit statistics...\n";
	node* curr_node = stats->headNode();
	char* buf_loc = 0;
	double* val = 0;
	size_t i = 0;
	size_t block = 1;
	for(countmap::iterator datIter = data->begin(); datIter != data->end(); ++datIter)
	{
		for (j = 0; j < set->size(); ++j)
		{
			buf_loc = curr_node->buf + sizeof(double) * i;
			val = reinterpret_cast<double*> (buf_loc);
			*val = calcWGOF(p[j], datIter->second.count, &(*set)[j]);
			++i;
			if (i >= stats->blockSize())
			{
				++block;
				if (block <= stats->nBlocks())
				{
					curr_node = curr_node->next;
					buf_loc = curr_node->buf;
					i = 0;
				}
			}
		}
	}

	// deallocate space for probability vector
	for (j = 0; j < set->size(); ++j)
		delete [] p[j];
}


// clear and deallocates "stat" member
bool kmer::clearStat ()
{
	if (stat != 0)
	{
		for (unsigned int i = 0; i < statsize; ++i)
			delete [] stat[i];
		delete [] stat;
		stat = 0;
		statsize = 0;
		return true;
	}
	else
	{
		fprintf(stderr, "WARNING: member stat is already empty\n");
		return false;
	}
}

// printCounts prints kmers and counts
void kmer::printCounts (std::ofstream& os, const countmap* kmers) const
{
	if (kmers->size() < 1)
	{
		fprintf(stderr, "No elements in kmer hash in call to kmer::printCounts\n");
		fail = 1;
		return;
	}
	std::stringstream seqstream;
	for (countmap::const_iterator kIter = kmers->begin(); kIter != kmers->end(); ++kIter)
	{
		os << numtoseq(kIter->first.id, seqstream);
		for (unsigned int k = 0; k < kIter->second.count.size(); ++k)
			os << "\t" << std::setw(12) << std::right << kIter->second.count[k];
	}
	os << "\n";
}

size_t kmer::nkmers ()
{
	return kmertypes;
}
