/*
 * parseData.h
 */

#ifndef PARSEDATA_H_
#define PARSEDATA_H_

#include <string>
#include <vector>
#include <fstream>

bool getFILE (std::fstream &, const char*, const char*);
int fexists (const char*);
std::vector<std::string> split (const std::string&, char);

#endif /* PARSEDATA_H_ */
