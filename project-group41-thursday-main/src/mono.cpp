#include "dy4.h"
#include "filter.h"
#include "fourier.h"
#include "genfunc.h"
#include "iofunc.h"
#include "logfunc.h"
#include <cmath>
#include <cstdlib>
#include <list>
#include <iostream>
#include <chrono>
#include <thread>
#include <queue>
#include <mutex>
#include <condition_variable>
void mono(std::queue<std::vector<float> > &my_queue, std::mutex &my_mutex, std::condition_variable &my_cvar, int mode){
    int rf_Fs = 2.4e6;
	int rf_Fc = 100e3;
	int rf_taps = 50;
	int rf_decim = 10;

    int rf1_Fs = 2.5e6;
	int mode0_Fs = 240e3;
	int mode0_decim = 5;
	int mode0_Fc =16e3;
	int mode0_taps = 50;

	int expand_rate = 24;
	int mode1_Fs = 250e3*expand_rate;
	int mode1_Fc = mode0_Fc;
	int mode1_decim = 125;
	int mode1_taps = 1250;

    std::vector<float> rf_coeff;
	std::vector<float> rf1_coeff;
	std::vector<float> mode0_coeff;
	std::vector<float> mode1_coeff;

	impulseResponseLPF(rf_Fs, rf_Fc, rf_taps, rf_coeff);
	impulseResponseLPF(rf1_Fs, rf_Fc, rf_taps, rf1_coeff);
	impulseResponseLPF(mode0_Fs, mode0_Fc, mode0_taps,mode0_coeff);
	impulseResponseLPF(mode1_Fs, mode1_Fc, mode1_taps,mode1_coeff);
	int BLOCK_SIZE = 1024*rf_decim*mode0_decim*2;
	std::vector<float> state_i_lpf_100k;
	state_i_lpf_100k.resize(rf_taps-1,0.0);
	std::vector<float> state_q_lpf_100k;
	state_q_lpf_100k.resize(rf_taps-1,0.0);
	std::vector<float> state_phase;
	state_phase.resize(2,0.0);
	std::vector<float> mode0_lpf_16k;
	std::vector<float> mode1_lpf_16k;
	mode0_lpf_16k.resize(mode0_taps-1,0.0);
	mode1_lpf_16k.resize(mode1_taps-1,0.0);
	std::vector<float> temp_I;
	temp_I.resize(BLOCK_SIZE/2,0.0);
	std::vector<float> temp_Q;
	temp_Q.resize(BLOCK_SIZE/2,0.0);
	std::vector<float> i_filt;
  std::vector<float> q_filt;
	std::vector<float> mode0_filt;
	std::vector<float> mode1_filt;
	std::vector<float> i_ds;
  std::vector<float> q_ds;
	std::vector<float> mode0_block;
	mode0_block.resize(BLOCK_SIZE/rf_decim/mode0_decim,0.0);
	std::vector<float> mode1_block;
	mode1_block.resize(BLOCK_SIZE*13/rf_decim/mode1_decim,0.0);
	std::vector<float> fm_demod;
	std::vector<float> exp_fm_demod;
	int mode0_data_size;
	std::vector<short int> mode0_data;
	mode0_data.resize(BLOCK_SIZE,0.0);
	int mode1_data_size;
	std::vector<short int> mode1_data;
	mode1_data.resize(BLOCK_SIZE,0.0);
	std::vector<float> block_data(BLOCK_SIZE);
	
	
	for (unsigned int block_id=0;; block_id++)
	{
		std::unique_lock<std::mutex> my_lock(my_mutex);
		if (my_queue.empty()){
			my_cvar.wait(my_lock);
		}
		std::vector<float> fm_demod = my_queue.front();
		my_queue.pop();
		std::cerr << "POPPED" << std::endl;
		if (mode==0)
		{
		 //convolveFIRblock(step, y, x, h, zi )
		 convolveFIRblock(mode0_decim, mode0_filt, fm_demod, mode0_coeff, mode0_lpf_16k);
		 downsample(mode0_filt, mode0_block, mode0_decim);
		 for (unsigned int n=0; n< mode0_block.size();n++){

			 if (std::isnan(mode0_block[n])){
				 mode0_data[n] = 0;
				 //mode0_data[(block_id*mode0_block.size())+n] = 0;
			 }
			 else{
				 //mode0_data[(block_id*mode0_block.size())+n] = static_cast<short int>(mode0_block[n]*16384);
				 mode0_data[n] = static_cast<short int>(mode0_block[n]*16384);

			 }
		 }
		 fwrite(&mode0_data[0],sizeof(short int),mode0_block.size(), stdout);
	 	}
		else if(mode==1)
		{
			expandBlock(exp_fm_demod, fm_demod, expand_rate);
			convolveFIRblock(mode1_decim, mode1_filt, exp_fm_demod, mode1_coeff, mode1_lpf_16k);
			downsample(mode1_filt, mode1_block, mode1_decim);
			for (unsigned int n=0; n< mode1_block.size();n++){

				if (std::isnan(mode1_block[n])){
					mode1_data[n] = 0;
				}
				else{
					mode1_data[n] = static_cast<short int>(mode1_block[n]*16384);
				}
			}
			fwrite(&mode1_data[0],sizeof(short int),mode1_block.size(), stdout);
		}
		my_lock.unlock();
		my_cvar.notify_one();
}
}
