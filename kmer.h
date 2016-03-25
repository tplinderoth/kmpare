/*
 * kmer.h
 */

#ifndef KMER_H_
#define KMER_H_

#include <vector>
#include <fstream>
#include <string>
#include <unordered_map>
#include <iostream>
#include <streambuf>
#include <sstream>
#include <iomanip>
#include "memPool.h"

template <class T>
class Array
{
public:
		T& operator[] (size_t i)
        {
				if (i >= size())
				{
					fprintf(stderr, "ERROR: subscript %lu out of range", i);
					exit(1);
				}
                return data[i];
        }

		T operator[] (size_t i) const
        {
				if (i >= size())
				{
					fprintf(stderr, "ERROR: subscript %lu out of range", i);
					exit(1);
				}
                return data[i];
        }

        void setSize(size_t size)
        {
				sz = size;
                data = new T[size];
                for(unsigned int long i = 0; i < size; ++i)
                	data[i] = 0;
        }

         size_t size() const
        {
        	 return sz;
        }

         Array ()
			 : data(0),
			   sz(0)
         { }

         Array ( const Array& oldarr)
			 : sz( oldarr.sz )
         {
        	 data = new T[sz];
        	 for( size_t i = 0; i < sz; ++i)
        		 data[i]= oldarr.data[i];
         }

        ~Array ()
        {
        	delete [] data;
        	data = 0;
        }


private:
        T* data;
        size_t sz;
};

template <class T>
struct Value
{
	Array<unsigned int> count;
	T* data;
};

template <class T>
struct Key
{
	Array<T> id;
	size_t* tableSize;

	bool operator==(const Key<T>& a) const
	{
		for(size_t i = 0; i < id.size(); ++i)
		{
			if (id[i] != a.id[i])
				return false;
		}
		return true;
	}
};

template <class T>
struct KeyHasher
{
	size_t operator() (const Key<T>& key) const
	{
		const int r = 31;
		unsigned long int h = 0;
		for(size_t i = 0; i < key.id.size(); ++i)
			h = (r * h + key.id[i]) % *(key.tableSize);
		return h;
	}
};

typedef std::unordered_map< Key<long int>, Value<double>, KeyHasher<long int> > countmap;

class kmer
{
public:
	//public functions
	kmer ();
	~kmer ();
	bool clearStat ();
	void parseJellyCounts (std::vector<std::string>& files);
	int jellyMerLength (std::ifstream& is);
	unsigned int long estLines (std::ifstream& is, int merlength, const int nonseq_n);
	std::string numtoseq (const Array<long int>& num, std::stringstream& ss) const;
	void fit(countmap* data, MemPool<double>* stats, std::vector< std::vector<unsigned int> >* set);
	template <class T> double calcWGOF (double p [], const Array<T>& obs, std::vector<unsigned int>* idx);
	template <class T> T arraySum (const Array<T>& v, std::vector<unsigned int>* index);
	void printCounts (std::ofstream& os, const countmap* kmers) const;
	template <class T> void printStats (std::ofstream& os, const countmap* kmers, size_t nstats, const MemPool<T>* stats) const;
	size_t nkmers ();
	// public data members
	mutable int fail;
	double** stat;
	size_t statsize;
	countmap datamap; // kmer-specific library counts
	Array<size_t> libtotal; // library-specific total counts across all kmers
	size_t kmerN;
private:
	//private functions
	void seqtonum (const std::string& s, Array<long int>& num, const int maxdigit);
	template <class T> int ndigit (T number);
	void libProbs (double p [], std::vector<unsigned int>* idx, Array<size_t>& lib_count);
	// private data members
	const int nonseq_char; // number of characters in each jellyfish file line, excluding the kmer, for estimating file size
	const float xtra_reserve; // allocates #_lines_in_1st_file * xtra_reserve more space for member "counts"
	unsigned int nlibs; // number of libraries to analyze
	size_t kmertypes; // number of actual different kmers in dataset
	size_t storage; // number of potential different kmer types to accommodate
};

// printStats prints kmers, counts, and additional information
template <class T> void kmer::printStats (std::ofstream& os, const countmap* kmers, size_t nstats, const MemPool<T>* stats) const
{
	if (kmers->size() < 1)
	{
		fprintf(stderr, "No elements in kmer hash in call to kmer::printStats\n");
		fail = 1;
		return;
	}
	if (stats->initial() == 0)
	{
		fprintf(stderr, "No statistics stored in memory pool in call to kmer::printStats\n");
		fail = 1;
		return;
	}

	std::stringstream seqstream;
	node* curr_node = stats->headNode();
	char* buf_loc = 0;
	size_t i = 0;
	size_t j = 0;
	size_t block = 1;
	T* val = 0;

	for (countmap::const_iterator kIter = kmers->begin(); kIter != kmers->end(); ++kIter)
	{
		os << numtoseq(kIter->first.id, seqstream);
		for (unsigned int k = 0; k < kIter->second.count.size(); ++k)
		{
			os << "\t" << std::setw(12) << std::right << kIter->second.count[k];
		}
		for (j = 0; j < nstats; ++j)
		{
			buf_loc = curr_node->buf + sizeof(T) * i;
			val = reinterpret_cast<T*> (buf_loc);
			os << "\t" << std::setw(12) << std::setprecision(5) << std::scientific << std::right << *val;
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
		os << "\n";
	}
}

#endif /* KMER_H_ */
