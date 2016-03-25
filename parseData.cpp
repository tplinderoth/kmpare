/*
 * parseData.cpp
 */

#include "parseData.h"
#include <sys/stat.h>
#include <sstream>
#include <cstring>

// getFILE is a wrapper for getting files
bool getFILE(std::fstream &fp, const char* fname, const char* mode)
{
	int writeFile = 0;
	if (strcmp(mode, "out") == 0)
	{
		writeFile = 1;

		if(writeFile && fexists(fname))
		{
			fprintf(stderr,"File already exists: %s\n",fname);
			return false;
		}

		fp.open(fname, std::ios::out);
	}
	else if (strcmp(mode, "app") == 0)
		fp.open(fname, std::ios::app);
	else if (strcmp(mode, "in") == 0)
		fp.open(fname, std::ios::in);

	if( !fp )
	{
		fprintf(stderr,"Error opening FILE handle for file: %s\n",fname);
		fp.close();
		return false;
	}
	return true;
}

// fexists finds out if a file exists
int fexists(const char* str)
{
	struct stat buffer;
 	return (stat(str, &buffer )==0 );
}

// splits a string based on a delimiter
std::vector<std::string> split (const std::string& s, char delim)
{
        std::vector<std::string> elems;
        std::stringstream ss(s);
        std::string sholder;
        while (std::getline(ss, sholder, delim))
        {
                if (!sholder.empty())
                        elems.push_back(sholder);
        }
        return elems;
}
