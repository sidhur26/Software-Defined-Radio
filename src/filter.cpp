/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#include "dy4.h"
#include "filter.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <complex>
#include <cmath>
#include <chrono>

#define PI 3.14159265358979323846

void delay_block(const std::vector<float> &input, std::vector<float> &output, std::vector<float> &state) {
    
    
    output.clear();
    state.clear();
    
    for (unsigned int i = 0; i < state.size(); i++) { 
        output.push_back(state[i]);
    }

    for (unsigned int i = 0; i < input.size() - state.size(); i++) { 
        output.push_back(input[i]);
    }

    state.assign(input.end() - state.size(), input.end());
}


void all_pass_filter(std::vector<float> &output, 
                        const std::vector<float> &input, 
                        std::vector<float> &state, 
                        int n_taps) {

    // Determine the number of previous state samples to retain for processing
    int stateSize = (n_taps - 1) / 2;
    
    // Prepare the output vector to hold the processed data
    output.clear();

    // Add the retained state samples to the output to maintain continuity
    output.insert(output.end(), state.begin(), state.end());

    // Append the new input samples to the output while excluding the last few
    // that are equal to the size of the state; these are stored for the next block
    output.insert(output.end(), input.begin(), input.end() - stateSize);

    // Update the state with the last samples from the current input block
    // These samples will be the state for the next block
    state.clear();
    state.insert(state.end(), input.end() - stateSize, input.end());
}



void fmPLL(std::vector<float> &input, float frequency, float sampling_frequency, float nco_scale, 
           float phase_adjust, float norm_bandwidth, std::vector<float> &PLL_state, std::vector<float> &nco_out) {

    // Constants for the phase-locked loop
    const float Cp = 2.666;
    const float Ci = 3.555;

    // Loop filter gains
    float Kp = norm_bandwidth * Cp;
    float Ki = norm_bandwidth * norm_bandwidth * Ci;
    
    // Ensure nco_out is large enough to hold all the outputs plus one for state continuity
    nco_out.resize(input.size() + 1, 0.0);

    // Extract state variables from the state vector
    float integrator = PLL_state[0];
    float phase_estimate = PLL_state[1];
    float feedback_I = PLL_state[2];
    float feedback_Q = PLL_state[3];
    float trig_offset = PLL_state[4];
    
    // The first element of nco_out is set to the last value of the previous block for continuity
    nco_out[0] = PLL_state[5];

    // Process each sample in the input signal
    for (size_t k = 0; k < input.size(); k++) {
        // Calculate errors by multiplying the input signal with the feedback
        float error_I = input[k] * feedback_I;
        float error_Q = input[k] * (-feedback_Q);

        // Calculate the phase detector error
        float error_D = std::atan2(error_Q, error_I);

        // Update the integrator and phase estimate
        integrator += Ki * error_D;
        phase_estimate += Kp * error_D + integrator;

        // Increment the trigonometric offset for the next sample
        trig_offset += 1;

        // Calculate the argument for trigonometric functions
        float trig_arg = 2 * M_PI * (frequency / sampling_frequency) * trig_offset + phase_estimate;
        
        // Update feedback components based on new phase
        feedback_I = std::cos(trig_arg);
        feedback_Q = std::sin(trig_arg);
        
        // Compute the NCO output for the current sample
        nco_out[k + 1] = std::cos(trig_arg * nco_scale + phase_adjust);
    }

    // Save current state back to the state vector for next block processing
    PLL_state[0] = integrator;
    PLL_state[1] = phase_estimate;
    PLL_state[2] = feedback_I;
    PLL_state[3] = feedback_Q;
    PLL_state[4] = trig_offset;
    PLL_state[5] = nco_out.back(); // Save the last NCO output value for continuity
}


void impulseResponseLPF(float Fs, float Fc, unsigned short int num_taps, std::vector<float> &h)
{
	// allocate memory for the impulse response
	h.clear(); h.resize(num_taps, 0.0);
    float cutoff = Fc / (Fs / 2);
 

    for (int i = 0; i < num_taps; i++) {
        if (i == (num_taps - 1) / 2) {
            h[i] = cutoff;
        } else {
                h[i] = cutoff * std::sin(PI * cutoff * (i - (num_taps - 1) / 2)) / (PI * cutoff * (i - (num_taps - 1) / 2));
            
        }
        h[i] *= std::pow(std::sin(i * PI / num_taps), 2);
    }
    
    
}

// function to compute the filtered output "y" by doing the convolution
// of the input data "x" with the impulse response "h"

void convolveFIR(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h)
{
    // allocate memory for the output (filtered) data
    y.clear();
    y.resize(x.size() + h.size() - 1, 0.0);

    const double epsilon = 1e-6;

    for (int n = 0; n < (int)y.size(); n++)
    {
        for (int k = 0; k < (int)h.size(); k++)
        {
            int idx = n - k;
            if (idx >= 0 && idx < (int)x.size())
            {
                y[n] += h[k] * x[idx];
            }
        }
        if (std::abs(y[n]) < epsilon)
        {
            y[n] = 0.0;
        }
    }
}


void convolveFIR_state(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h, std::vector<float> &state) {
    
    y.clear();
    y.resize(x.size(), 0.0);
    unsigned int size_state = h.size() - 1;

    // Performing convolution
    for (int n = 0; n < (int)x.size(); n++) {
        for (int k = 0; k < (int)h.size(); k++) {
        if (n-k >= 0) {
            y[n] += h[k] * x[n - k]; }
        else { 
            y[n] += h[k] * state[size_state + n - k]; 
            
            }
        }
    }

    // Updating the state for the next block
    state.clear();
    state.insert(state.end(), x.end() - size_state, x.end()); 
}

void impulseResponseBPF(const float Fs,
                        const float Fb, 
                        const float Fe, 
                        const unsigned short int n_taps, 
                        std::vector<float> &filter_coefficients)
{
    // Clearing the existing vector and resizing it to the new number of taps, aand initialize with zeros
    filter_coefficients.clear();
    filter_coefficients.resize(n_taps, 0.0);

    //float numerator_denominator; 

    // Calculating normalized center frequency and bandwidth of the filter
    float normCenter = ((Fe + Fb) / 2) / (Fs / 2);
    float normPass = (Fe - Fb) / (Fs / 2);

    // Loop through each filter tap to calculate its coefficient
    for (int i = 0; i < n_taps - 1; i++){
        if (i == (n_taps - 1) / 2){
            // Assign normalized passband width to the center coefficient
            filter_coefficients[i] = normPass;
        }
        else{
            // Calculate the sinc function for non-center coefficients
            //numerator_denominator = (PI * (normPass / 2) * (i - (n_taps - 1) / 2));
            // Avoid division by zero, using sinc function
            filter_coefficients[i] = normPass * (std::sin((PI * (normPass / 2) * (i - (n_taps - 1) / 2))) / (PI * (normPass / 2) * (i - (n_taps - 1) / 2)));
        }
        // Shift the filter to the correct center frequency and apply window
        filter_coefficients[i] *= std::cos((i - (n_taps - 1) / 2) * PI * normCenter); // Shift center frequency
        //filter_coefficients[i] *= std::cos(i  * PI * normCenter); // Shift center frequency
    
        filter_coefficients[i] *= std::pow(std::sin((i * PI) / n_taps), 2); // Apply window function
    }
}

//Function converts to raw data normalizes
void readStdinBlockData(unsigned int num_samples,
                        unsigned int block_id,
                        std::vector<float> &block_data)
{
    std::vector<char> raw_data(num_samples);
    std::cin.read(reinterpret_cast<char*>(&raw_data[0]), num_samples*sizeof(char));
    for (int k=0; k<(int)num_samples; k++) {
        // automatically normalizes the data to the range -1 to +1
        block_data[k] = float(((unsigned char)raw_data[k]-128)/128.0);
    }
}


void decim2(std::vector<float> &output, 
            const std::vector<float> &input, 
            const std::vector<float> &filter_coefficients, 
            std::vector<float> &state, 
            int downsample) {

    output.clear(); 
    output.resize(input.size() / downsample, 0.0); // Resize output vector to appropriate length based on downsampling rate
	
	
    for (int i = 0; i < (int)input.size() / downsample; i++) {
		//int down_index = i * downsample;
        for (int j = 0; j < (int) filter_coefficients.size(); j++) {
            // Downsample the signal
            if(i * downsample - j >= 0){ 
                // Perform the convolution using filter coefficients
                output[i] += filter_coefficients[j] * input[i * downsample - j];
            }
            else {
                // Use the saved state for indices that are negative after downsampling
                output[i] += filter_coefficients[j] * state[(int)state.size() + (i * downsample - j)];
            }
        }
    }

    // Update the state for the next block of processing
    state.assign(input.end() - filter_coefficients.size() + 1, input.end());
}


void decim4(std::vector<float> &output, 
            const std::vector<float> &input, 
            const std::vector<float> &filter_coefficients, 
            std::vector<float> &state, 
            int downsample, 
            int upsample) {


	int input_size = input.size();
    int state_vector_size = state.size();
    
    // Resize the output vector to account for both upsampling and downsampling
    output.clear();
    output.resize(input_size * upsample / downsample, 0.0);
   

    // Iterate over the resized output
    for (int i = 0; i < (int)output.size(); i++) {
        // Calculate the index in the input corresponding to the current output element
        int input_index = (i * downsample - upsample) / upsample;

        for (int j = 0; j < (int) filter_coefficients.size(); j++) {
            // Calculate the actual index in the input or state to be used for convolution
            int index = input_index + j;

            if(index >= 0 && index < input_size) {
                // Perform the convolution using filter coefficients on input
                output[i] += filter_coefficients[j] * input[index];
            } else if (index < 0) {
                // Use the saved state for indices that are negative
                int state_index = state.size() + index;
                if (state_index >= 0 && state_index < state_vector_size) {
                    output[i] += filter_coefficients[j] * state[state_index];
                }
            }
        }
    }

    // Update the state for the next block of processing
    int state_size = std::min((int)input.size(), (int)filter_coefficients.size() - 1);
    state.assign(input.end() - state_size, input.end());
}


void fmDemodulation(const std::vector<float> &I, const std::vector<float> &Q, 
                    float &prev_I, float &prev_Q, std::vector<float> &fm_demod) {
    // Allocating memory for output
    fm_demod.clear();
    fm_demod.resize(I.size());

    for (unsigned int k = 0; k < I.size(); k++) {
        if (k == 0) {
            fm_demod[k] = 0.0; // Initializing first element
        } else {
            // FM demodulation formula
            fm_demod[k] = (I[k] * (Q[k] - prev_Q) - Q[k] * (I[k] - prev_I)) / ((I[k] * I[k]) + (Q[k] * Q[k]));
        }

        // State saving vars
        prev_I = I[k];
        prev_Q = Q[k];
    }
}

void stereo_combiner(std::vector<float> &stereo_output, const std::vector<float> &mono_input, const std::vector<float> &stereo_input) {
	
	stereo_output.clear();
	int input_size = mono_input.size();
	stereo_output.resize(input_size*2);
	
	for (unsigned int j = 0; j < mono_input.size(); j++){
		stereo_output[2*j] = (mono_input[j] + stereo_input[j])/2.0;
		stereo_output[2*j+1] = (mono_input[j] - stereo_input[j])/2.0;
	}

}

void stereo_mixer(std::vector<float> &mixed_stereo, const std::vector<float> &stereo_input, const std::vector<float> &carrier_input) {
	
	int input_size = stereo_input.size();
	mixed_stereo.resize(input_size);
	
	for (int i = 0; i < input_size; i++) {
		mixed_stereo[i] = stereo_input[i] * carrier_input[i] *2;
	}

}

void squaring_nonlinearity(const std::vector<float> &input, std::vector<float> &squared) {
    // Ensure the squared signal vector has the same size as the input signal
    squared.resize(input.size());

    // Perform the squaring operation on each element of the input signal
    for (size_t i = 0; i < input.size(); i++) {
        squared[i] = input[i] * input[i];
    }
}

