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
#include <chrono>
#include <mutex>
#include <condition_variable>
#include <atomic>


// declaration of a function prototypes
void delay_block(const std::vector<float> &inblock, std::vector<float> &outblock, std::vector<float> &stateblock);
void all_pass_filter(std::vector<float> &output, const std::vector<float> &input, std::vector<float> &state, int n_taps);
void fmPLL(std::vector<float> &input, float frequency, float sampling_frequency, float nco_scale, float phase_adjust, float norm_bandwidth, std::vector<float> &PLL_state, std::vector<float> &nco_out);
void impulseResponseLPF(float, float, unsigned short int, std::vector<float> &);
void convolveFIR(std::vector<float> &, const std::vector<float> &, const std::vector<float> &);
void convolveFIR_state(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h, std::vector<float> &state);
void convolutionFIR_state(const std::vector<float>& audioData, const std::vector<float>& filterCoeff, const std::vector<float>& zi, std::vector<float>& filteredSig, std::vector<float>& newState);
void impulseResponseBPF(const float Fs,const float Fb, const float Fe, const unsigned short int n_taps, std::vector<float> &filter_coefficients);
void readStdinBlockData(unsigned int num_samples, unsigned int block_id, std::vector<float> &block_data);
void decim2(std::vector<float> &output, const std::vector<float> &input, const std::vector<float> &filter_coefficients, std::vector<float> &state, int downsample);
void decim4(std::vector<float> &output, const std::vector<float> &input, const std::vector<float> &filter_coefficients, std::vector<float> &state, int downsample, int upsample);
void fmDemodulation(const std::vector<float> &I, const std::vector<float> &Q, float &prev_I, float &prev_Q, std::vector<float> &fm_demod);
void stereo_combiner(std::vector<float> &stereo_output, const std::vector<float> &mono_input, const std::vector<float> &stereo_input);
void stereo_mixer(std::vector<float> &mixed_stereo, const std::vector<float> &stereo_input, const std::vector<float> &carrier_input);
void squaring_nonlinearity(const std::vector<float> &input, std::vector<float> &squared);
                        
                        
                        
                        
                        


#endif // DY4_FILTER_H
