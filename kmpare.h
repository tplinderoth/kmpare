/*
 * kmpare.h
 */

#ifndef KMPARE_H_
#define KMPARE_H_

#include <vector>
#include <fstream>

// version
const char * version = "0.1.1"; // 7 December 2014

// functions
bool parseArgs (int argc, char** argv, std::vector<std::string>* ifname, std::vector< std::vector<unsigned int> >* cmpindex, std::string& ofname);
std::vector<unsigned int> parseSet (int argc, char** argv, int& pos);
void printHeader (std::ofstream& os, unsigned int nlibs, const std::vector< std::vector<unsigned int> >* sets);
void info (const char* v);

#endif /* KMPARE_H_ */
