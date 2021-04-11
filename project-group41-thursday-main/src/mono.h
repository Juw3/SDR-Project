/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#ifndef DY4_MONO_H
#define DY4_MONO_H

// add headers as needed
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <complex>
#include <queue>
#include <mutex>
#include <condition_variable>
#include <thread>
// declaration of a function prototypes
void mono(std::queue<std::vector<float> > &my_queue, std::mutex &my_mutex, std::condition_variable &my_cvar, int mode);
#endif // DY4_MONO_H
