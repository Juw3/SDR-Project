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
#include "rds.h"
#include "mono.h"
#include "stereo.h"
#include "rf_frontend.h"


int main(int argc, char* argv[]){
	std::queue<std::vector <float> > my_queue;
	std::mutex my_mutex;
	std::condition_variable my_cvar;
	// std::queue<std::vector <float> > my_queue_rds;
	// std::mutex my_mutex_rds;
	// std::condition_variable my_cvar_rds;
	
	int mode = 0;
	// int stereo_done = 0;
	// int rds_done = 0;
	// int stereo_index = 0;
	// int rds_index = 0;
	if(argc < 2){
		std::cerr << " Operating in default mode 0" << std::endl;
	} else if (argc == 2){
		mode = atoi(argv[1]);
		if (mode != 1){
			std::cerr << "Wrong mode" << mode << std::endl;
			exit(1);
		}
	}else{
		std::cerr << "Usage: " << argv[0] << std::endl;
		std::cerr << "or" <<std::endl;
		std::cerr << "Usage: " << argv[0] << " 1" << std::endl;
		exit(1);
	}
	
	// std::thread t_mono  = std::thread(mono,std::ref(my_queue),std::ref(my_mutex),std::ref(my_cvar), mode);
	std::thread t_rf = std::thread(rf_frontend, std::ref(my_queue), std::ref(my_mutex), std::ref(my_cvar), mode);

	std::thread t_mono_stereo  = std::thread(mono_stereo,std::ref(my_queue),std::ref(my_mutex),std::ref(my_cvar), mode);

	// std::thread t_rf = std::thread(rf_frontend, std::ref(my_queue), std::ref(my_mutex), std::ref(my_cvar), std::ref(my_queue_rds),std::ref(my_mutex_rds),std::ref(my_cvar_rds), mode);

	// t_rds = std::thread(rds,std::ref(my_queue_rds), std::ref(my_mutex_rds), std::ref(my_cvar_rds));

	// if (mode == 0){
	// 	t_rf = std::thread(rf_frontend, std::ref(my_queue), std::ref(my_mutex), std::ref(my_cvar), std::ref(my_queue_rds),std::ref(my_mutex_rds),std::ref(my_cvar_rds), mode);

	// 	t_rds = std::thread(rds,std::ref(my_queue_rds), std::ref(my_mutex_rds), std::ref(my_cvar_rds));
	// }
	
	
	t_rf.join();
	// t_mono.join();
	t_mono_stereo.join();

	// if (mode == 0){
	// 	t_rds.join();
	// }

    // mono(mode);
 return 0;

}

/////////////////////////////TO RUN RDS/////////////////////////////////
////////////////COMMENT LINES ABOVE//////////////////////////////
////////////////UNCOMMENT LINES BELOW/////////////////
// #include "dy4.h"
// #include "filter.h"
// #include "fourier.h"
// #include "genfunc.h"
// #include "iofunc.h"
// #include "logfunc.h"
// #include <cmath>
// #include <cstdlib>
// #include <list>
// #include <iostream>
// #include <chrono>
// #include <thread>
// #include <fstream>
// #include <vector>
// #include <complex>

// int main(){
//   float rf_Fs = 2.4e6;
//   float rf_Fc = 100e3;
//   int rf_taps = 151;
//   int rf_decim = 10;
//   float ncoScale = 0.5;
//   float phaseAdjust = 0;
//   float normBandwidth = 0.01;
//   float feedbackI = 1;
//   float feedbackQ = 0;
//   float integrator = 0;
//   float phaseEst = 0;
//   float last_I = 1;
// 	float last_Q = 0;
//   float trigOffset = 0;
//   float fb0 = 54e3;
//   float fe0 = 60e3;
//   float fb1 = 113.5e3;
//   float fe1 = 114.5e3;
//   float fs_stereo = 240e3;
//   float freq = 114e3;
// 	int BLOCK_SIZE = 1024*50;
//   int mini = 0;
//   int minq = 0;
//   int posi = 0;
//   int posq = 0;
//   int mark_i=0;
//   int mark_q=0;
// 	std::vector<float> rf_coeff;
// 	std::vector<float> rds0_coeff;
// 	std::vector<float> rds1_coeff;
// 	std::vector<float> rds1_l_coeff;
// 	std::vector<float> rds2_l_coeff;
// 	std::vector<float> impulseRRC;
// 	std::vector<float> state_i_lpf_100k;
// 	state_i_lpf_100k.resize(rf_taps-1,0.0);
// 	std::vector<float> state_q_lpf_100k;
// 	state_q_lpf_100k.resize(rf_taps-1,0.0);
// 	std::vector<float> state_phase;
// 	state_phase.resize(2,0.0);
// 	std::vector<float> temp_I;
// 	temp_I.resize(BLOCK_SIZE/2,0.0);
// 	std::vector<float> temp_Q;
// 	temp_Q.resize(BLOCK_SIZE/2,0.0);
// 	std::vector<float> i_filt;
//   std::vector<float> q_filt;
// 	std::vector<float> i_ds;
//   std::vector<float> q_ds;
// 	std::vector<float> fm_demod;
// 	std::vector<float> rds0_bpf_16k;
// 	rds0_bpf_16k.resize(150,0.0);
// 	std::vector<float> rds1_bpf_16k;
//   rds1_bpf_16k.resize(150,0.0);
// 	std::vector<float> rds0_i_lpf_16k;
//   rds0_i_lpf_16k.resize(150,0.0);
// 	std::vector<float> rds0_q_lpf_16k;
//   rds0_q_lpf_16k.resize(150,0.0);
// 	std::vector<float> rds1_i_lpf_16k;
//   rds1_i_lpf_16k.resize(19*151-1,0.0);
// 	std::vector<float> rds1_q_lpf_16k;
//   rds1_q_lpf_16k.resize(19*151-1,0.0);
// 	std::vector<float> rds_i_rrc_16k;
//   rds_i_rrc_16k.resize(150,0.0);
// 	std::vector<float> rds_q_rrc_16k;
//   rds_q_rrc_16k.resize(150,0.0);
//   impulseResponseLPF(2.4e6, 100e3, 151, rf_coeff);
// 	bandPass(54e3, 60e3, 240e3, 151, rds0_coeff);
//   bandPass(113.5e3, 114.5e3, 240e3, 151, rds1_coeff);
// 	impulseResponseLPF(240e3, 3e3, 151, rds1_l_coeff);
//   impulseResponseLPF(4.56e6, 19e3, 2869, rds2_l_coeff);
//   impulseResponseRootRaisedCosine(57000, 151, impulseRRC);
// 	std::vector<float> rds0_filt;
// 	int rds_decim = 1;
// 	std::vector<float> nl_sqr;
// 	std::vector<float> rds1_filt;
// 	std::vector<float> I_Out;
// 	std::vector<float> Q_Out;
// 	std::vector<float> mixer_I;
// 	std::vector<float> mixer_Q;
// 	std::vector<float> rds0_l_i_filt;
// 	std::vector<float> rds0_l_q_filt;
// 	std::vector<float> rds_i_expand;
// 	std::vector<float> rds_q_expand;
//   std::vector<float> rds1_l_i_filt;
// 	std::vector<float> rds1_l_q_filt;
//   std::vector<float> rds_ds_i;
//   std::vector<float> rds_ds_q;
//   std::vector<float> rds_rrc_i_filt;
//   std::vector<float> rds_rrc_q_filt;
//   std::vector<float> rrc_rec_I;
//   std::vector<float> rrc_rec_Q;
//   std::vector<std::vector<int> > impulsevec{{1000000000},{0100000000},{0010000000},{0001000000},{0000100000},{0000010000},{0000001000},{0000000100},{0000000010},{0000000001},{1011011100},{0101101110},{0010110111},{1010000111},{1110011111},{1100010011},{1101010101},{1101110110},{0110111011},{1000000001},{1111011100},{0111101110},{0011110111},{1010100111},{1110001111},{1100011011}};
//   int lastvalue = 0;
//   std::vector<int> ddcode;
//   //std::vecotr<int> getddcode;
//   int matrixsum;
//   // std::vector<float> saving_state_i;
//   // std::vector<float> saving_state_q;
// 	std::vector<float> block_data(BLOCK_SIZE);
//   std::vector<int> savestate;
// 	auto start_time = std::chrono::high_resolution_clock::now();
// 	for (unsigned int block_id=0;; block_id++)
// 	{

// 		readStdinBlockData(BLOCK_SIZE, block_id, block_data);
// 		if((std::cin.rdstate())!=0){
// 				 std::cerr << "End of input stream reached" <<std::endl;
// 			 //	for (auto j=0;j<BLOCK_SIZE;j+=10){std::cerr << "block data" << block_data[j] << std::endl;}
// 				 //std::cerr << "Printing audio block size" << mode0_block.size() << std::endl;
// 				// fwrite(&mode0_stereo_data[0],sizeof(short int),mode0_stereo_data.size(), stdout);
// 				 auto stop_time = std::chrono::high_resolution_clock::now();
// 				 std::chrono::duration<double, std::milli> run_time = stop_time - start_time;
// 				 std::cerr << "effBlockConv ran for " << run_time.count() << " milliseconds" << std::endl;
// 				 exit(1);
// 		 }
// 	 		// std::cerr << "Read block" << block_id << std::endl;
// 	 		int j = 0;
// 	 		for(int i = 0; i<BLOCK_SIZE; i=i+2)
// 	 		{
// 	 				temp_I[j] = block_data[i];
// 	 				temp_Q[j] = block_data[i+1];
// 	 				j++;
// 	 	  }
// 			convolveFIRblock(rf_decim, i_filt, temp_I, rf_coeff, state_i_lpf_100k);
// 			convolveFIRblock(rf_decim, q_filt, temp_Q, rf_coeff, state_q_lpf_100k);
// 			//downsample(y,x,ratio)
// 			//y is input, x is output of downsample
// 			downsample(i_filt, i_ds, rf_decim);
// 			downsample(q_filt, q_ds, rf_decim);
// 			//fmDemodArctan(I, Q, prev_phase,fmdemod)
// 			fmDemodArctan(i_ds, q_ds, state_phase, fm_demod);
//       //bpf with Fb = 54KHz, Fe = 60KHz, Fs = 240KHz
// 		 	convolveFIRblock(rds_decim, rds0_filt, fm_demod, rds0_coeff, rds0_bpf_16k);
// 	    //squaring
// 			nl_sqr.resize(rds0_filt.size(), 0.0);
// 			for(auto i = 0; i<rds0_filt.size(); i++)
// 			{
// 				 nl_sqr[i] = rds0_filt[i]*rds0_filt[i];
// 			}
//       //bpf with Fb = 113.5KHz, Fe = 114.5KHz, Fs = 240KHz
// 			convolveFIRblock(rds_decim, rds1_filt, nl_sqr, rds1_coeff, rds1_bpf_16k);
//       //fmpll
// 			fmPllRDS(rds1_filt, 114e3, 240e3, ncoScale, phaseAdjust, normBandwidth, I_Out, Q_Out, feedbackI, feedbackQ, last_I, last_Q, integrator, phaseEst, trigOffset);
//       //mixer
// 			mixer_I.resize(rds0_filt.size(), 0.0);
// 			mixer_Q.resize(rds0_filt.size(), 0.0);
// 			for(auto k = 0; k<rds0_filt.size(); k++)
// 			{
// 				mixer_I[k] = I_Out[k]*rds0_filt[k];
// 				mixer_Q[k] = Q_Out[k]*rds0_filt[k];
// 			}
//       //lpf with fc=3KHz fs= 240KHz
// 			convolveFIRblock(rds_decim, rds0_l_i_filt, mixer_I, rds1_l_coeff, rds0_i_lpf_16k);
// 			convolveFIRblock(rds_decim, rds0_l_q_filt, mixer_Q, rds1_l_coeff, rds0_q_lpf_16k);
//       //upsampler with ratio 19

// 			expandBlock(rds_i_expand, rds0_l_i_filt, 19);
// 			expandBlock(rds_q_expand, rds0_l_q_filt, 19);
//       //lpf with Fc=19KHZ, Fs=19*240KHz, Ntaps=19*151
// 			convolveFIRblock(80, rds1_l_i_filt, rds_i_expand, rds2_l_coeff, rds1_i_lpf_16k);
// 			convolveFIRblock(80, rds1_l_q_filt, rds_q_expand, rds2_l_coeff, rds1_q_lpf_16k);
//       //downsampler with ratio 80
//       downsample(rds1_l_i_filt, rds_ds_i, 80);
//       downsample(rds1_l_q_filt, rds_ds_q, 80);
//       //rrc

//       convolveFIRblock(rds_decim, rds_rrc_i_filt, rds_ds_i, impulseRRC, rds_i_rrc_16k);
// 	  convolveFIRblock(rds_decim, rds_rrc_q_filt, rds_ds_q, impulseRRC, rds_q_rrc_16k);

//       // for(auto i =0; i<rds_rrc_i_filt.size();i++)
//       // {
//       //   std::cerr << "Printing I" << rds_rrc_i_filt[i] << std::endl;
//       //   std::cerr << "Printing Q " <<rds_rrc_q_filt[i] << std::endl;
//       // }
//       std::vector<float> t1;
//       t1.resize(rds_rrc_i_filt.size(),0.0);
//       for(auto i=0; i<t1.size(); i++)
//       {
//         t1[i] = i;
//       }
//       std::vector<float> t0;
//       t0.resize(rds_ds_i.size(),0.0);
//       for(auto i=0; i<t0.size(); i++)
//       {
//         t0[i] = i;
//       }
//       // int max_index_I = 0;
//       // int max_index_Q = 0;
//       // for(auto i = 0; i< 24; i++){
//       //   if (abs(rds_rrc_i_filt[max_index_I])<abs(rds_rrc_i_filt[i]))
//       //       max_index_I = i;
//       //   if (abs(rds_rrc_q_filt[max_index_Q])<abs(rds_rrc_q_filt[i]))
//       //       max_index_Q = i;
//       // }
//       // std::vector<float> rrc_rec_I;
//       // std::vector<float> rrc_rec_Q;
//       // std::cerr << "Printing I" << posi << std::endl;
//       // std::cerr << "Printing Q " <<posq << std::endl;
//       // rrc_rec_I.resize(rds_rrc_i_filt.size()/24,0.0);
//       // rrc_rec_Q.resize(rds_rrc_q_filt.size()/24,0.0);
//       // int mini = 0;
//       // int minq = 0;
//       // int posi = 0;
//       // // int posq = 0;
//       // int posi = 0;
//       // int posq = 0;
//       // for(auto i =0; i<200; i++)
//       //  {
//       //    if (rds_rrc_i_filt[posi]<rds_rrc_i_filt[i])
//       //    {
//       //      posi = i;
//       //    }
//       //    if (rds_rrc_q_filt[posq]<rds_rrc_q_filt[i])
//       //    {
//       //      posq = i;
//       //    }
//       //  }
//       //int rem = rds_rrc_i_filt.size()%12;
//       if (block_id==0)
//       {
//       for(auto i =0; i<24; i++)
//       {
//         if (rds_rrc_i_filt[mini]>rds_rrc_i_filt[i])
//         {
//           mini = i;
//         }
//       }
//       if (mini > 23)
//       {
//         posi = 12;
//       }
//       else if(mini <= 23)
//       {
//         posi = 0;
//       }
//       posq = 12;
//       std::cerr << "Printing I" << posi << std::endl;
//       std::cerr << "Printing Q " << posq << std::endl;
//      }
//      // int rem = rds_rrc_i_filt.size()%24;
//      // posi = posi+rem;
//      // posq = posq+rem;
//      // std::vector<float> saving_state_i;
//      //saving_state_i.resize(rem,0.0);
//      // std::vector<float> saving_state_q;
//      //saving_state_q.resize(rem,0.0);
//      rrc_Recovery(posi,posq,mark_i,mark_q, rds_rrc_i_filt, rds_rrc_q_filt, rrc_rec_I, rrc_rec_Q);
//      std::vector<int> manchester;
//      manchester.resize(rrc_rec_I.size()/2,0.0);
//      int index = 0;
//      for (auto i = 0; i<rrc_rec_I.size();i=i+2)
//      {
//        if (rrc_rec_I[i]>0&&rrc_rec_I[i+1]<0)
//        {
//          manchester[index] = 1;
//          index++;
//        }
//        else if (rrc_rec_I[i]<0&&rrc_rec_I[i+1]>0)
//        {
//          manchester[index] = 0;
//          index++;
//        }
//      }
//      std::vector<int> ddcode;
//      int matrix_pos = 0;
//     ddecoding(manchester, lastvalue, ddcode);
//     // for (auto f=0; f<10; f++)
//     // {
//     //   std::cerr << "Printing frame" <<ddcode[f] << std::endl;
//     // }
//     if(block_id==0){
//     savestate.resize(ddcode.size()%26, 0.0);
//    }
//     // std::cout << "Printing size" <<ddcode.size() << std::endl;
//     std::vector<int> frame_out;
//     frame(ddcode,impulsevec,frame_out, savestate);

//     // for (auto f=0; f<10; f++)
//     // {
//     //   std::cerr << "Printing frame" <<frame_out[f] << std::endl;
//     // }

//        // std::cerr << "Printing I" <<mark_i << std::endl;
//        // std::cerr << "Printing Q " <<mark_q << std::endl;
//     // matrixsum=0;
//     // for(auto i=0;i < BLOCK_SIZE;){
//     //   for(auto k=0;k<26;k++){
//     //       for(auto j=0;j<10;j++){
//     //         matrixsum+=ddcode[i]*implusevec[k][j];
//     //       }
//     //       getddcode[k]=matrixsum;
//     //       i++;
//     //   }
//     // }

//     // if (block_id == 0){
//     //
//      if (block_id == 0){
//        for (auto f=0; f<10; f++)
//        {
//          std::cerr << "Printing frame" <<frame_out[f] << std::endl;
//        }
//         logVector("rrc", rrc_rec_I, rrc_rec_Q); // log only positive freq
//         logVector("demod_time_O", t1, rds_ds_i);
//         logVector("demod_time_I", t0, rds_rrc_i_filt);
//         logVector("demod_time_O1", t0, rds_rrc_q_filt);
//         std::cerr << "Run: gnuplot -e 'set terminal  png size 1024,768' example.gnuplot > ../data/example.png\n";
//     //}
// }
// }
// }
