/*
Comp Eng 3DY4 (Computer Systems Integration Project)
Copyright by Nicola Nicolici
Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

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
void mono_stereo(std::queue<std::vector<float> > &my_queue, std::mutex &my_mutex, std::condition_variable &my_cvar, int mode)
{

	int rf_Fs = 2.4e6;
	int rf_Fc = 100e3;
	int rf_taps = 101;
	int rf_decim = 10;

    int rf1_Fs = 2.5e6;
	int mode0_Fs = 240e3;
	int mode0_decim = 5;
	int mode0_Fc =16e3;
	int mode0_taps = 101;
    int BLOCK_SIZE = 1024*rf_decim*mode0_decim*2;
	int expand_rate = 24;
	int mode1_Fs = 250e3*expand_rate;
	int mode1_Fc = mode0_Fc;
	int mode1_decim = 125;
	int mode1_taps = 2424;

	float fb0 = 18.5e3;
	float fe0 = 19.5e3;
	float fb1 = 22e3;
	float fe1 = 54e3;
	float fs0_stereo = 240e3;
	float fs1_stereo = 250e3;
	float ncoScale = 2;
	float phaseAdjust = 0;
	float normBandwidth = 0.01;
	int stereo_taps = 101;
	float freq = 19e3;
	float feedbackI=1;
	float feedbackQ=0;
	float lastnco=1;
	float integrator=0;
	float phaseEst=0;
	float trigOffset=0;
	int stereo_decim = 1;

	std::vector<float> rf_coeff;
	std::vector<float> rf1_coeff;
	std::vector<float> mode0_coeff;
	std::vector<float> mode1_coeff;
	std::vector<float> stereo00_coeff;
	std::vector<float> stereo01_coeff;
	std::vector<float> stereo10_coeff;
	std::vector<float> stereo11_coeff;
	std::vector<float> pllin;
	std::vector<float> ncoOut;
	std::vector<float> mixer;
	std::vector<float> stereo0_lpf_16k;
	std::vector<float> stereo1_lpf_16k;
	stereo0_lpf_16k.resize(stereo_taps-1,0.0);
	stereo1_lpf_16k.resize(stereo_taps-1,0.0);
	std::vector<float> stereo11_lpf_16k;
	stereo11_lpf_16k.resize(2424-1,0.0);
	std::vector<float> stereo0_filt;
	std::vector<float> stereo1_filt;
	std::vector<float> exp_ncoOut;
	std::vector<float> ncoOutfilt;
	std::vector<float> nco_lpf_16k;
	nco_lpf_16k.resize(mode1_taps-1,0.0);
	std::vector<float> nco1_block;
	nco1_block.resize(BLOCK_SIZE*12/rf_decim/mode1_decim,0.0);
	std::vector<float> nco0_block;
	nco0_block.resize(BLOCK_SIZE/2/rf_decim/mode0_decim,0.0);
	 std::vector<float> LAC;
	 std::vector<float> RAC;
	impulseResponseLPF(mode0_Fs, mode0_Fc, mode0_taps,mode0_coeff);
	impulseResponseLPF(mode1_Fs, mode1_Fc, mode1_taps,mode1_coeff);
	bandPass(fb0,fe0,fs0_stereo,stereo_taps,stereo00_coeff);
	bandPass(fb1,fe1,fs0_stereo,stereo_taps,stereo01_coeff);
	bandPass(fb0,fe0,fs1_stereo,stereo_taps,stereo10_coeff);
	bandPass(fb1,fe1,fs1_stereo,stereo_taps,stereo11_coeff);

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
	mode0_block.resize(BLOCK_SIZE/2/rf_decim/mode0_decim,0.0);
	std::vector<float> mode1_block;
	mode1_block.resize(BLOCK_SIZE*12/rf_decim/mode1_decim,0.0);
	std::vector<float> fm_demod;
	std::vector<float> exp_fm_demod;
	std::vector<float> mode0_data;
	mode0_data.resize(926720,0.0);
	std::vector<float> mode1_data;
	mode1_data.resize(821788,0.0);
	std::vector<short int> mode0_stereo_data;
	mode0_stereo_data.resize(2*926720,0.0);
	std::vector<short int> mode1_stereo_data;
	mode1_stereo_data.resize(2*821788,0.0);
	std::vector<float> mode0_stereo_block;
	mode0_stereo_block.resize(BLOCK_SIZE/mode0_decim/rf_decim,0.0);
	std::vector<float> mode1_stereo_block;
	mode1_stereo_block.resize(BLOCK_SIZE*24/mode1_decim/rf_decim,0.0);
	std::vector<float> block_data(BLOCK_SIZE);
	auto start_time = std::chrono::high_resolution_clock::now();
	for (unsigned int block_id=0;; block_id++)
	{
        std::unique_lock<std::mutex> my_lock(my_mutex);
		if (my_queue.empty()){
			my_cvar.wait(my_lock);
		}
		std::vector<float> fm_demod = my_queue.front();
		my_queue.pop();
		if (mode==0)
		 {

		 convolveFIRblock(mode0_decim, mode0_filt, fm_demod, mode0_coeff, mode0_lpf_16k);
		 downsample(mode0_filt, mode0_block, mode0_decim);
		 		convolveFIRblock(stereo_decim, stereo0_filt, fm_demod, stereo00_coeff, stereo0_lpf_16k);
		 		convolveFIRblock(stereo_decim, stereo1_filt, fm_demod, stereo01_coeff, stereo1_lpf_16k);
		 		fmPll(stereo0_filt, freq, fs0_stereo, ncoScale, phaseAdjust, normBandwidth, ncoOut, feedbackI, feedbackQ, lastnco, integrator, phaseEst, trigOffset);
		 	  mixer.resize(stereo0_filt.size(),0.0);
		 			for(auto k=0; k<stereo0_filt.size(); k++)
		 			{
		 				mixer[k] = 2*ncoOut[k]*stereo1_filt[k];
		 			}
		 			convolveFIRblock(mode0_decim, ncoOutfilt, mixer, mode0_coeff, mode0_lpf_16k);
		 		  downsample(ncoOutfilt, nco0_block, mode0_decim);
		 			LAC.resize(nco0_block.size(),0.0);
		 			RAC.resize(nco0_block.size(),0.0);
		 			//audio_block.resize(2*nco_block.size(),0.0);
		 			for(auto z=0; z<nco0_block.size(); z++)
		 			{
		 				LAC[z]= (nco0_block[z]+mode0_block[z])/2;
		 				RAC[z]= (mode0_block[z]-nco0_block[z])/2;
		 			}
		 			int f = 0;
		 		  for(auto x=0; x<mode0_stereo_block.size(); x=x+2)
		 			{
		 				mode0_stereo_block[x] = LAC[f];
		 				mode0_stereo_block[x+1] = RAC[f];
		 				f=f+1;
		 			}
		 		for (unsigned int n=0; n< mode0_stereo_block.size();n++){

		 			if (std::isnan(mode0_stereo_block[n])){
		 				mode0_stereo_data[(block_id*mode0_stereo_block.size())+n] = 0;
		 			}
		 			else{
		 				mode0_stereo_data[(block_id*mode0_stereo_block.size())+n] = static_cast<short int>(mode0_stereo_block[n]*16384);
		 			}
		 		}
	 //
		  fwrite(&mode0_stereo_data[block_id*mode0_stereo_block.size()],sizeof(short int),mode0_stereo_block.size(), stdout);
}
		else if(mode==1)
		{
			expandBlock(exp_fm_demod, fm_demod, expand_rate);
			//convolveFIRblock(mode1_decim, mode1_filt, exp_fm_demod, mode1_coeff, mode1_lpf_16k);
			convolveFIRblock_mode1(mode1_decim,expand_rate ,mode1_filt, exp_fm_demod, mode1_coeff, mode1_lpf_16k);
			downsample(mode1_filt, mode1_block, mode1_decim);
			convolveFIRblock(stereo_decim, stereo0_filt, fm_demod, stereo10_coeff, stereo0_lpf_16k);
			convolveFIRblock(stereo_decim, stereo1_filt, fm_demod, stereo11_coeff, stereo1_lpf_16k);
			fmPll(stereo0_filt, freq, fs1_stereo, ncoScale, phaseAdjust, normBandwidth, ncoOut, feedbackI, feedbackQ, lastnco, integrator, phaseEst, trigOffset);
			mixer.resize(stereo0_filt.size(),0.0);
			for(auto k=0; k<stereo0_filt.size(); k++)
			{
				mixer[k] = 2*ncoOut[k] * stereo1_filt[k];
			}
			expandBlock(exp_ncoOut, mixer, expand_rate);
			//convolveFIRblock(mode1_decim, ncoOutfilt, exp_ncoOut, mode1_coeff, stereo11_lpf_16k);
			convolveFIRblock_mode1(mode1_decim, expand_rate, ncoOutfilt, exp_ncoOut, mode1_coeff, stereo11_lpf_16k);
			downsample(ncoOutfilt, nco1_block, mode1_decim);
			LAC.resize(nco1_block.size(),0.0);
			RAC.resize(nco1_block.size(),0.0);
			//audio_block.resize(2*nco_block.size(),0.0);
			for(auto z=0; z<nco1_block.size(); z++)
			{
				LAC[z]= (nco1_block[z]+mode1_block[z])/2;
				RAC[z]= (mode1_block[z]-nco1_block[z])/2;
			}
			int f = 0;
			for(auto x=0; x<mode1_stereo_block.size(); x=x+2)
			{
				mode1_stereo_block[x] = LAC[f];
				mode1_stereo_block[x+1] = RAC[f];
				f=f+1;
			}
			for (unsigned int n=0; n< mode1_stereo_block.size();n++){

				if (std::isnan(mode1_stereo_block[n])){
					mode1_stereo_data[(block_id*mode1_stereo_block.size())+n] = 0;
				}
				else{
					mode1_stereo_data[(block_id*mode1_stereo_block.size())+n] = static_cast<short int>(mode1_stereo_block[n]*16384);
				}
			}
 //
		fwrite(&mode1_stereo_data[block_id*mode1_stereo_block.size()],sizeof(short int),mode1_stereo_block.size(), stdout);

	    }
        std::cerr << "Write block" << block_id << std::endl;
        // if (block_id == 1024){
        //     stop_playing = 0;
        // }

        my_lock.unlock();
		my_cvar.notify_one();
    }
}
