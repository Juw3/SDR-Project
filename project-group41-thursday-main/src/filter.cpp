/*
Comp Eng 3DY4 (Computer Systems Integration Project)

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
#include <math.h>

void fmPll(std::vector<float> &pllIn, float freq, float Fs, float ncoScale, float phaseAdjust, float normBandwidth, std::vector<float> &ncoOut, float &feedbackI, float &feedbackQ, float &lastnco, float &integrator, float &phaseEst, float &trigOffset)
{
	float Cp = 2.666;
	float Ci = 3.555;

	float Kp = (normBandwidth)*Cp;
	float Ki = (normBandwidth * normBandwidth) * Ci;

	ncoOut.resize(pllIn.size() + 1, 0.0);
	ncoOut[0] = lastnco;
	float errorI;
	float errorQ;
	float errorD;
	float trigArg;

	int k = 0;
	for (k = 0; k < pllIn.size(); k++)
	{
		errorI = pllIn[k] * feedbackI;
		errorQ = pllIn[k] * (-1) * feedbackQ;

		errorD = atan2(errorQ, errorI);

		integrator = integrator + Ki * errorD;

		phaseEst = phaseEst + Kp * errorD + integrator;

		trigArg = 2 * PI * (freq / Fs) * (trigOffset + k + 1) + phaseEst;
		feedbackI = cos(trigArg);
		feedbackQ = sin(trigArg);
		lastnco = cos(trigArg * ncoScale + phaseAdjust);
		ncoOut[k + 1] = lastnco;
	}
	trigOffset += pllIn.size();
}

void bandPass(float f_b, float f_e, float f_s, unsigned short int N_taps, std::vector<float> &h)
{
	h.resize(N_taps, 0.0);
	float norm_Center = (((f_e + f_b) / 2) / (f_s / 2));
	float norm_Pass = (f_e - f_b) / (f_s / 2);
	unsigned short int N = int(N_taps);
	for (auto i = 0; i < N; i++)
	{
		if (i == (N - 1) / 2)
		{
			h[i] = norm_Pass;
		}
		else
		{
			h[i] = norm_Pass * (std::sin(PI * (norm_Pass / 2) * (i - (N - 1) / 2)) / (PI * (norm_Pass / 2) * (i - (N - 1) / 2)));
		}
		h[i] *= std::cos(i * PI * norm_Center);
		h[i] *= std::sin(i * PI / N) * std::sin(i * PI / N);
	}
}
// function to compute the impulse response "h" based on the sinc function
void expandBlock(std::vector<float> &yb, std::vector<float> &xb, int step)
{
	yb.resize(step * xb.size(), 0.0);
	int n = 0;
	for (auto i = 0; i < xb.size(); i++)
	{
		yb[n] = xb[i] * 24;
		n = n + step;
	}
}
// void expandBlock(std::vector<float> &yb, std::vector<float> &xb, int step)
// {
//     yb.resize(step * xb.size(), 0.0);
//     int n = 0;
//     for (auto i = 0; i < xb.size(); i++)
//     {
//         yb[n] = 24 * xb[i];
//         n = n + step;
//     }
// }
// int queue_block[QUEUE_BLOCKS][BLOCK_SIZE];
// template <size_t size_x, size_t size_y>
void fmPllRDS(std::vector<float> &pllIn, float freq, float Fs, float ncoScale, float phaseAdjust, float normBandwidth, std::vector<float> &I_Out, std::vector<float> &Q_Out, float &feedbackI, float &feedbackQ, float &last_I, float &last_Q, float &integrator, float &phaseEst, float &trigOffset)
{
    float Cp = 2.666;
    float Ci = 3.555;

    float Kp = (normBandwidth)*Cp;
    float Ki = (normBandwidth * normBandwidth) * Ci;

    I_Out.resize(pllIn.size() + 1, 0.0);
    Q_Out.resize(pllIn.size() + 1, 0.0);
    I_Out[0] = last_I;
    Q_Out[0] = last_Q;
    float errorI;
    float errorQ;
    float errorD;
    float trigArg;

    int k = 0;
    for (k = 0; k < pllIn.size(); k++)
    {
        errorI = pllIn[k] * feedbackI;
        errorQ = pllIn[k] * (-1) * feedbackQ;

        errorD = atan2(errorQ, errorI);

        integrator = integrator + Ki * errorD;

        phaseEst = phaseEst + Kp * errorD + integrator;

        trigArg = 2 * PI * (freq / Fs) * (trigOffset + k + 1) + phaseEst;
        feedbackI = cos(trigArg);
        feedbackQ = sin(trigArg);
        last_I = cos(trigArg * ncoScale + phaseAdjust);
        last_Q = sin(trigArg * ncoScale + phaseAdjust);
        I_Out[k + 1] = last_I;
        Q_Out[k + 1] = last_Q;
    }
    trigOffset += pllIn.size();
}
void impulseResponseRootRaisedCosine(float Fs, int N_taps, std::vector<float> &impulseResponseRRC)
{

    //duation for each symbol - do NOT be changed for RDS!
    float T_symbol = 1 / 2375.0;

    //roll-off factor (greater than 0 and smaller than 1)
    float beta = 0.90;
    float value = T_symbol / (4 * beta);
    float t = 0;
    //the RRC inpulse response that will be computed in this function
    impulseResponseRRC.resize(N_taps, 0.0);

    for (auto k = 0; k < N_taps; k++)
    {
        t = (k - N_taps / 2) / Fs;
        //we ignore the 1/T_symbol scale factor
        if (t == 0.0)
        {
            impulseResponseRRC[k] = 1.0 + beta * ((4 / PI) - 1);
        }
        else if (abs(t) == value)
        {
            impulseResponseRRC[k] = (beta / sqrt(2)) * (((1 + 2 / PI) * (sin(PI / (4 * beta)))) + ((1 - 2 / PI) * (cos(PI / (4 * beta)))));
        }
        else
        {
            impulseResponseRRC[k] = (sin(PI * t * (1 - beta) / T_symbol) + 4 * beta * (t / T_symbol) * cos(PI * t * (1 + beta) / T_symbol)) / (PI * t * (1 - (4 * beta * t / T_symbol) * (4 * beta * t / T_symbol)) / T_symbol);
        }
    }

}
void frame(std::vector<int> &ddcode, std::vector<std::vector<int> > &impulsevec, std::vector<int> &frame_out, std::vector<int> &savestate)
{
    std::vector<int> frame;
    frame.resize(26, 0);
    frame_out.resize(10, 0);
    int len = savestate.size();
    for (auto i = 0; i < len; i++)
    {
        frame[i] = savestate[i];
    }
    if (len < 26)
    {
        for (auto j = 0; j < 26 - len; j++)
        {
            frame[j + len] = ddcode[j];
        }
        savestate.resize(26 - len, 0.0);
        for (auto j = 0; j < 26 - len; j++)
        {
            savestate[j] = ddcode[j + len];
        }
        for (auto k = 0; k < 10; k++)
        {
            for (auto z = 0; z < 26; z++)
            {
                frame_out[k] += frame[z] * impulsevec[z][k];
            }
        }
    }
    else
    {
        savestate.resize(0, 0);
        for (auto k = 0; k < 10; k++)
        {
            for (auto z = 0; z < 26; z++)
            {
                frame_out[k] += frame[z] * impulsevec[z][k];
            }
        }
    }
}
void ddecoding(std::vector<int> manchester, int &lastvalue, std::vector<int> &ddcode)
{
    std::vector<int> manch_buf;
    manch_buf.resize(manchester.size() + 1, 0.0);
    manch_buf[0] = lastvalue;
    lastvalue = manchester[manchester.size() - 1];
    ddcode.resize(manchester.size(), 0.0);
    for (auto i = 0; i < manch_buf.size() - 1; i++)
    {
        manch_buf[i + 1] = manchester[i];
    }
    for (auto i = 0; i < ddcode.size(); i++)
    {
        ddcode[i] = ((!manch_buf[i]) * manch_buf[i + 1]) + (manch_buf[i] * (!manch_buf[i + 1]));
    }
}
void rrc_Recovery(int &posi, int &posq, int &mark_i, int &mark_q, std::vector<float> &rds_rrc_i_filt, std::vector<float> &rds_rrc_q_filt, std::vector<float> &rrc_rec_I, std::vector<float> &rrc_rec_Q)
{
    int tempi = posi;
    int tempq = posq;
    rrc_rec_I.resize(rds_rrc_i_filt.size() / 24, 0.0);
    rrc_rec_Q.resize(rds_rrc_q_filt.size() / 24, 0.0);
    for (auto i = 0; i < rrc_rec_I.size(); i++)
    {
        rrc_rec_I[i] = rds_rrc_i_filt[tempi];
        rrc_rec_Q[i] = rds_rrc_q_filt[tempq];
        tempi = tempi + 24;
        tempq = tempq + 24;
        // std::cerr << "Printing I" << i << std::endl;
        // std::cerr << "Printing Q " << maxq << std::endl;
    }
    mark_i = (rds_rrc_i_filt.size() - posi) % 24;
    mark_q = (rds_rrc_q_filt.size() - posq) % 24;
    posi = 24 - mark_i;
    posq = 24 - mark_q;
}
void fmDemodArctan(std::vector<float> &I, std::vector<float> &Q, std::vector<float> &prev_phase, std::vector<float> &fmdemod)
{

	fmdemod.resize(I.size(), 0.0);
	std::vector<float> Iext, Qext;
	Iext.resize(I.size() + 1, 0.0);
	Qext.resize(Q.size() + 1, 0.0);
	Iext[0] = prev_phase[0];
	Qext[0] = prev_phase[1];
	for (auto j = 1; j < Iext.size(); j++)
	{
		Iext[j] = I[j - 1];
		Qext[j] = Q[j - 1];
	}
	int k;

	for (k = 0; k < Iext.size() - 1; k++)
	{
		float sqI = pow(I[k], 2);
		float sqQ = pow(Q[k], 2);
		if (sqI + sqQ == 0)
		{
			fmdemod[k] = 0;
			// queue_block[queue_entry][k] = 0;
		}
		else
		{
			fmdemod[k] = (Iext[k + 1] * (Qext[k + 1] - Qext[k]) - Qext[k + 1] * (Iext[k + 1] - Iext[k])) / (sqI + sqQ);
			// queue_block[queue_entry][k] = (Iext[k+1]*(Qext[k+1]-Qext[k])-Qext[k+1]*(Iext[k+1]-Iext[k]))/(sqI+sqQ);
		}
	}
	prev_phase[0] = Iext[k];
	prev_phase[1] = Qext[k];
}
void estimatePSD(std::vector<float> &bin_data, int Fs, std::vector<float> &freq, std::vector<float> &psd_est)
{
	int freq_bins = NFFT;
	float df = Fs / freq_bins;
	freq.resize(int(freq_bins / 2), 0.0);
	for (auto i = 0; i < (Fs / 2) / df; i++)
	{
		freq[i] = i * df;
		freq[i] = freq[i] / 1000;
	}

	std::vector<float> hann;
	hann.resize(freq_bins, 0.0);
	for (auto j = 0; j < hann.size(); j++)
	{
		hann[j] = pow(sin(j * PI / freq_bins), 2);
	}
	std::vector<float> psd_list;
	int no_segments = floor(bin_data.size() / freq_bins);
	psd_list.resize(freq_bins / 2 + no_segments * freq_bins / 2, 0.0);
	for (auto k = 0; k < no_segments; k++)
	{
		//apply the hann window (using pointwise multiplication)
		//before computing the Fourier transform on a segment
		std::vector<float> temp;
		temp.assign(bin_data.begin() + k * freq_bins, bin_data.begin() + (k + 1) * freq_bins);
		std::vector<float> windowed_samples;
		windowed_samples.resize(freq_bins, 0.0);
		for (auto i = 0; i < windowed_samples.size(); i++)
		{
			windowed_samples[i] = temp[i] * hann[i];
		}
		std::vector<std::complex<float> > Yf;
		//std::vector <std::complex<float>> Yf;
		DFT(windowed_samples, Yf);

		Yf.assign(Yf.begin(), Yf.begin() + int(freq_bins / 2));
		std::vector<float> psd_seg;
		psd_seg.resize(freq_bins / 2, 0.0);
		for (auto i = 0; i < psd_seg.size(); i++)
		{
			psd_seg[i] = 0.00001628 * pow(abs(Yf[i]), 2);
			psd_seg[i] = 2 * psd_seg[i];
			psd_seg[i] = 10 * log10(psd_seg[i]);
		}
		psd_list.insert(psd_list.begin() + k * psd_seg.size(), psd_seg.begin(), psd_seg.end());
	}
	psd_est.resize(int(freq_bins / 2), 0.0);
	for (auto n = 0; n < int(freq_bins / 2); n++)
	{
		for (auto l = 0; l < no_segments; l++)
		{
			psd_est[n] += psd_list[n + l * int(freq_bins / 2)];
		}
		psd_est[n] = psd_est[n] / no_segments;
	}
}
void downsample(std::vector<float> &y, std::vector<float> &x, int ratio)
{
	x.resize(y.size() / ratio, 0.0);
	int j = 0;
	for (auto i = 0; i < y.size(); i = i + ratio)
	{
		x[j] = y[i];
		j++;
	}
}
void impulseResponseLPF(float Fs, float Fc, unsigned short int num_taps, std::vector<float> &h)
{
	// bring your own functionality
	// allocate memory for the impulse response
	h.resize(num_taps, 0.0);
	float norm = (Fc / (Fs / 2));
	for (auto i = 0; i < num_taps; i++)
	{
		if (i == ((num_taps - 1) / 2))
		{
			h[i] = norm;
		}
		else
		{
			h[i] = norm * (sin(PI * norm * (i - (num_taps - 1) / 2))) / (PI * norm * (i - (num_taps - 1) / 2));
		}
		h[i] = h[i] * sin(i * PI / num_taps) * sin(i * PI / num_taps);
	}
}
void convolveFIRblock(int step, std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h, std::vector<float> &zi)
{
	y.resize(x.size(), 0.0);
	// std::cout << "x.size: " << x.size() << " \n";
	// std::cout << "zi.size: " << zi.size() << " \n";
	for (int i = 0; i < x.size(); i = i + step)
	{ // loop through the size of x
		float sum = 0;
		for (int j = 0; j < h.size(); j++)
		{ //then loop through the size of h
			if ((i - j) >= 0)
			{ //if the indexes are both positive
				sum += h[j] * x[i - j];
			}
			else
			{ //if not then use values from zi
				// std::cout << "index: " << zi.size()+(i-j) << "\n";
				sum += h[j] * zi[zi.size() + (i - j)];
			}
		}
		y[i] = sum; //add value to y index
	}
	//calculating the new zi from x
	auto first = x.cbegin() + (x.size() - h.size() + 1);
	auto last = x.cbegin() + (x.size() - 1);
	std::vector<float> newZi(first, last);
	zi = newZi;
}

void convolveFIRblock_mode1(int step,int expand,std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h,std::vector<float> &zi)
{
  y.resize(x.size(), 0.0);

  for (auto i = 0; i < x.size(); i=i+step){ //every 'steps' block is downsampled and needs to be read
		float sum = 0;
    int start = (i/expand)*expand; //determine the start value

		for(auto j = start; i-j< h.size(); j-=expand){
      if(j >= 0){
        sum += x[j]*h[i-j]; //x contains mostly 0's from alot of 0s from upsampling, so only want to take every 24 block
      }
      else{ //if not then use values from zi
      	sum += h[i-j]*zi[zi.size()+ j];
      }
    }
		y[i] = sum;
  }

	//calculating the new zi from x
	auto first = x.cbegin() + (x.size()-h.size()+1);
	auto last = x.cbegin() + (x.size()-1);
	std::vector<float> newZi(first, last);
	zi = newZi;

}

// function to compute the filtered output "y" by doing the convolution
// of the input data "x" with the impulse response "h"
