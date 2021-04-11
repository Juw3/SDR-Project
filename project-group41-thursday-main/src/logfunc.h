/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#ifndef DY4_LOGFUNC_H
#define DY4_LOGFUNC_H

// add headers as needed
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <complex>

// declaration of a function prototypes


void logVector(const std::string filename, const std::vector<float> &x, const std::vector<float> &y);

void genIndexVector(std::vector<float> &x, const int size);

#endif // DY4_LOGFUNC_H
