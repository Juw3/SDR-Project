#ifndef DY4_RF_FRONTEND_H
#define DY4_RF_FRONTEND_H

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

void rf_frontend(std::queue<std::vector<float> > &my_queue, std::mutex &my_mutex, std::condition_variable &my_cvar, int mode);
#endif // DY4_RF_FRONTEND_H

