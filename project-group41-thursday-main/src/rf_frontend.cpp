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

void rf_frontend(std::queue<std::vector<float> > &my_queue, std::mutex &my_mutex, std::condition_variable &my_cvar, int mode){
    int rf_Fs = 2.4e6;
	int rf_Fc = 100e3;
	int rf_taps = 50;
	int rf_decim = 10;

	int rf1_Fs = 2.5e6;
	int mode0_decim = 5;
	int mode0_taps = 50;

	int expand_rate = 24;
	int mode1_Fs = 250e3 * expand_rate;
	int mode1_decim = 125;
	int mode1_taps = 1250;

	std::vector<float> rf_coeff, rf1_coeff, mode0_coeff, mode1_coeff;

	impulseResponseLPF(rf_Fs, rf_Fc, rf_taps, rf_coeff);
	impulseResponseLPF(rf1_Fs, rf_Fc, rf_taps, rf1_coeff);
	int BLOCK_SIZE = 1024 * rf_decim * mode0_decim * 2;
	std::vector<float> state_i_lpf_100k, state_q_lpf_100k;
	state_i_lpf_100k.resize(rf_taps - 1, 0.0);
	state_q_lpf_100k.resize(rf_taps - 1, 0.0);
	std::vector<float> state_phase;
	state_phase.resize(2, 0.0);
	std::vector<float> temp_I, temp_Q;
	temp_I.resize(BLOCK_SIZE / 2, 0.0);
	temp_Q.resize(BLOCK_SIZE / 2, 0.0);
	std::vector<float> i_filt, q_filt, mode0_filt, mode1_filt, i_ds, q_ds;
	std::vector<float> mode0_block, mode1_block;
	mode0_block.resize(BLOCK_SIZE /2/ rf_decim / mode0_decim, 0.0);
	mode1_block.resize(BLOCK_SIZE * 12 / rf_decim / mode1_decim, 0.0);
	std::vector<float> fm_demod;
	std::vector<float> block_data(BLOCK_SIZE);
    int QUEUE_BLOCKS = 5;
    int block_id = 0;
    int queue_block[QUEUE_BLOCKS][BLOCK_SIZE];
    for (unsigned int block_id = 0;; block_id++)
	{
        // unsigned int queue_entry = block_id % QUEUE_BLOCKS;
        readStdinBlockData(BLOCK_SIZE, block_id, block_data);
		if ((std::cin.rdstate()) != 0)
		{
            std::cerr << "End of input stream reached" << std::endl;
            //---------------------------------------------------------------------------------------------------
            //905 blocks:
            //writes to 312 BLOCK_SIZE
            //completes in 7000 seconds

            std::this_thread::sleep_for(std::chrono::milliseconds(21000));
            //	for (auto j=0;j<BLOCK_SIZE;j+=10){std::cerr << "block data" << block_data[j] << std::endl;}
            int mode_block = mode?mode1_block.size():mode0_block.size();
            std::cerr << "Printing audio block size" << mode_block << std::endl;
            //fwrite(&audio_data[0],sizeof(short int),audio_data.size(), stdout);
            // std::chrono::duration<double, std::milli> run_time = stop_time - start_time;
            // std::cerr << "effBlockConv ran for " << run_time.count() << " milliseconds" << std::endl;
            exit(1);
		}
        std::cerr << "Read block" << block_id << std::endl;
        int j = 0;
        for (int i = 0; i < BLOCK_SIZE; i = i + 2)
        {
            temp_I[j] = block_data[i];
            temp_Q[j] = block_data[i + 1];
            j++;
        }
        convolveFIRblock(rf_decim, i_filt, temp_I, rf_coeff, state_i_lpf_100k);
        convolveFIRblock(rf_decim, q_filt, temp_Q, rf_coeff, state_q_lpf_100k);
        downsample(i_filt, i_ds, rf_decim);
        downsample(q_filt, q_ds, rf_decim);
        fmDemodArctan(i_ds, q_ds, state_phase, fm_demod);

        std::unique_lock<std::mutex> my_lock(my_mutex);
        if (my_queue.size() == QUEUE_BLOCKS-1){
            my_cvar.wait(my_lock);
        }
        my_queue.push(fm_demod);
        // std::cerr << "Pushed!!" << std::endl;
        // my_queue.push((void *)&queue_block[queue_entry][0]);
        my_lock.unlock();
        my_cvar.notify_one();
        // block_id++;

        // std::unique_lock<std::mutex> my_lock_rds(my_mutex_rds);
        // if (my_queue_rds.size() == QUEUE_BLOCKS-1){
        //     my_cvar_rds.wait(my_lock);
        // }
        // my_queue_rds.push(fm_demod);
        // // std::cerr << "Pushed!!" << std::endl;
        // // my_queue.push((void *)&queue_block[queue_entry][0]);
        // my_lock_rds.unlock();
        // my_cvar_rds.notify_one();

    }
}
