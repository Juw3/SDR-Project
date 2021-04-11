/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#ifndef DY4_FILTER_H
#define DY4_FILTER_H

// add headers as needed
#include <iostream>
#include <vector>
#include <complex>
#include <list>

// declaration of a function prototypes
// template<size_t size_x, size_t size_y>
// size_t size_x, size_y;
void frame(std::vector<int> &ddcode, std::vector<std::vector<int> > &impulsevec, std::vector<int> &frame_out, std::vector<int> &savestate);
void ddecoding(std::vector<int> manchester, int &lastvalue, std::vector<int> &ddcode);
void rrc_Recovery(int &posi,int &posq, int &mark_i, int &mark_q, std::vector<float> &rds_rrc_i_filt, std::vector<float> &rds_rrc_q_filt, std::vector<float> &rrc_rec_I, std::vector<float> &rrc_rec_Q);
void fmPllRDS(std::vector<float> &pllIn, float freq, float Fs, float ncoScale, float phaseAdjust, float normBandwidth, std::vector<float> &I_Out, std::vector<float> &Q_Out,float &feedbackI, float &feedbackQ, float &last_I, float &last_Q, float &integrator, float &phaseEst, float &trigOffset);
void impulseResponseRootRaisedCosine(float Fs, int N_taps, std::vector<float> &impulseResponseRRC);

void fmPll(std::vector<float> &pllIn, float freq, float Fs, float ncoScale, float phaseAdjust, float normBandwidth, std::vector<float> &ncoOut, float &feedbackI, float &feedbackQ, float &lastnco, float &integrator, float &phaseEst, float &trigOffset);
void bandPass(float f_b, float f_e, float f_s, unsigned short int N_taps, std::vector<float> &h);
void expandBlock(std::vector<float> &yb, std::vector<float> &xb, int step);
void fmDemodArctan(std::vector<float> &I, std::vector<float> &Q, std::vector<float> &prev_phase, std::vector<float> &fmdemod);
void estimatePSD(std::vector<float> &bin_data, int Fs, std::vector<float> &freq,  std::vector<float> &psd_est);
void downsample(std::vector<float> &x, std::vector<float> &y, int ratio);
void impulseResponseLPF(float Fs, float Fc, unsigned short int num_taps, std::vector<float> &h);
void convolveFIRblock(int step, std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h, std::vector<float> &zi);
void convolveFIR(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h);
void convolveFIRblock_mode1(int step,int expand,std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h,std::vector<float> &zi);

#endif // DY4_FILTER_H
