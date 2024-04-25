
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


//GLOBAL VARIABLES
int BLOCK_SIZE;

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


void decim(std::vector<float> &output, const std::vector<float> &input, const int rf_decim){
	
    output.clear(); 
    for(unsigned int i=0; i<input.size(); i+=rf_decim){
           output.push_back(input[i]);
        
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
/*
void decim2(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h, std::vector<float> &filter_state, int down){ 
    y.clear(); y.resize(x.size()/down, 0.0); 

    for (int i = 0; i < (int)y.size(); i++) {
        for (int j = 0; j < (int) h.size(); j++) {
            //downsampling
            if(i*down-j >= 0){ 
                y[i] += h[j] * x[i*down-j];
            }
            //else if (-(i*down-j) <= filter_state.size()){ 
            else {
                y[i] += h[j] * filter_state[(int) filter_state.size() + (i*down - j)];
            }

        }
    }
    filter_state.assign(x.end() - (h.size()-1), x.end());
    //filter_state.assign(x.end() - h.size()+1, x.end());  
    
}
*/


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

int main(int argc, char* argv[])
{
	//------------------------------------------------------------MODE SELECTION---------------------------------------------------------------------
		// using namespace std::chrono;
		//default mode 0 
		int mode = 0;
		
		if (argc<2){
			std::cerr << "Operating in default mode 0" << std::endl;
			
		} else if (argc==2){
			mode=atoi(argv[1]);
			
			if (mode>3){
				std::cerr << "Wrong mode " << mode << std::endl;
				exit(1);
			}
		} else {
			std::cerr << "Usage: " << argv[0] << std::endl;
			std::cerr << "or " << std::endl;
			std::cerr << "Usage: " << argv[0] << " <mode>" << std::endl;
			std::cerr << "\t\t <mode> is a value from 0-3" << argv[0] << std::endl;
			exit(1);
		}

		std::cerr << "Operating in mode " << mode << std::endl;
	
        
        //Mode based variable declarations
        float rf_Fs;
        float rf_decim;
        float audio_decim;
        int audio_up;
        unsigned int BLOCK_SIZE;
        float audio_IF;
        
        if (mode == 0){
		audio_IF = 2.4e5;
        rf_Fs = 2.4e6;
        rf_decim = 10;
        audio_decim = 5;
        audio_up = 1;
        BLOCK_SIZE = 1024*rf_decim*((float)audio_decim/audio_up) * 2; 
        }
        
		else if (mode == 1){
		audio_IF = 1.2e5;
		rf_Fs = 9.6e5;
        rf_decim = 8;
        audio_decim = 3;
        audio_up = 1;
        BLOCK_SIZE = 1015*rf_decim*((float)audio_decim/audio_up) * 2; 
		}
		
		else if (mode == 2){
		audio_IF = 2.4e5;
		rf_Fs = 2.4e6;
        rf_decim = 10;
        audio_decim = 800;
        audio_up = 147;
        BLOCK_SIZE = 1030*rf_decim*((float)audio_decim/audio_up) * 2;
		}
		
		else {
		audio_IF = 3.6e5;
		rf_Fs = 1.44e6;
        rf_decim = 4;
        audio_decim = 400;
        audio_up = 49;
        BLOCK_SIZE = 1313*rf_decim*((float)audio_decim/audio_up) * 2;
		}
		
		// -------------------------------------------------- VARIABLE DECLARATIONS ---------------------------------------------------------
		// -------------------- MONO ------------------------
		float rf_Fc = 100e3;
		float rf_taps = 101;
		float audio_taps = 101;
		float audio_Fc = 16e3;
		
		std::vector<float>block_data(BLOCK_SIZE);
		std::vector<float>I(block_data.size() /2);
        std::vector<float>Q(block_data.size() /2);	
        std::vector<float> I_down;
        std::vector<float> Q_down;	
		float prev_I = 0.0;
        float prev_Q = 0.0;

		std::vector<float> rf_h;
		std::vector<float> audio_h;
        
        std::vector<float> IF; 
        std::vector<float> IF_filtered;
        std::vector<float> IF_down(BLOCK_SIZE);
        std::vector<float> IF_state(audio_taps - 1, 0.0);
        
		std::vector<float> I_state(rf_taps - 1, 0.0);
        std::vector<float> Q_state(rf_taps - 1, 0.0); 
        
        std::vector<float> I_filtered;
        std::vector<float> Q_filtered;
        
        //std::vector<float> mono_output;
       
        // -------------------- STEREO -----------------------
        /*
        float n_taps = 101;
        

        std::vector<float> extraction_h;
        std::vector<float> extraction_filtered;
        
        std::vector<float> carrier_h;
        std::vector<float> carrier_filtered;
        
        std::vector<float> carrier_state(rf_taps - 1, 0.0);
        std::vector<float> extraction_state(rf_taps - 1, 0.0);
        
        std::vector<float> delay_state(rf_taps - 1, 0.0);
        std::vector<float> IF_delayed;
        
        float PLL_frequency = 19e3;
        std::vector<float> mixed_stereo(BLOCK_SIZE/(2*rf_decim), 0.0);
        
        std::vector<float> stereo_state(n_taps, 0.0);
        std::vector<float> stereo_conv;
		std::vector<float> stereo_h;
		std::vector<float> stereo_output;
		
		std::vector<float> stereo_left_right(BLOCK_SIZE/(rf_decim*((float)audio_decim/audio_up)), 0.0);
		
		float integrator = 0.0;
		float phase_estimate = 0.0;
		float feedbackI = 1.0;
		float feedbackQ = 0.0; 
		float trig_offset = 0.0; 
		float nco_scale = 1; 
		float phase_adjust = 0.0; 
		float norm_bandwith = 0.001;
		std::vector<float> nco_out(BLOCK_SIZE/(2*rf_decim), 1.0);
		
		float extraction_fb = 22e3;
        float extraction_fe = 54e3;
        float stereo_decim = 1.0;
        float carrier_fb = 18.5e3;
        float carrier_fe = 19.5e3;
        
        //impulseResponseLPF(rf_Fs,rf_Fc,rf_taps,rf_h);
       // impulseResponseBPF(audio_IF, extraction_fb, extraction_fe, n_taps, extraction_h);
       // impulseResponseBPF(audio_IF, carrier_fb, carrier_fe, n_taps, carrier_h);
        //impulseResponseLPF(audio_IF,audio_Fc,rf_taps, audio_h);
       // impulseResponseLPF(audio_IF, audio_Fc, n_taps, stereo_h);
        
        
        std::cerr << "Size of rf_h: " << rf_h.size() << std::endl;
        std::cerr << "Size of extraction_h: " << extraction_h.size() << std::endl;
        std::cerr << "Size of carrier_h: " << carrier_h.size() << std::endl;
        std::cerr << "Size of audio_h: " << audio_h.size() << std::endl;
        std::cerr << "Size of stereo_h: " << stereo_h.size() << std::endl;
        */
        
        
        // ------------------ BLOCK READING ----------------------
        
       for(unsigned int block_id = 0; ; block_id++){
		
		//Read raw data 
		readStdinBlockData(BLOCK_SIZE, block_id, block_data);
		
		if((std::cin.rdstate()) != 0){
            std::cerr << "End of input stream reached!" << std::endl;
            exit(1);
        }
        
        //std::cerr << "Processing block " << block_id <<std::endl;  //output  message indicating that a block has been successfully read
        
		
		// ----------------------------------------------------------- READING -------------------------------------------------------------------------------
		//Splitting block data -> I, Q
		
        int iq_index = 0;
        for (unsigned int i=0; i < block_data.size(); i+=2) {
            I[iq_index] = block_data[i];
            Q[iq_index] = block_data[i + 1];
            iq_index++;
        }
        
        // -------------- RF FRONT END ------------------------------
        
        impulseResponseLPF(rf_Fs,rf_Fc,rf_taps,rf_h);
 
        //Decimation
        decim2(I_down, I, rf_h, I_state, rf_decim);        
        decim2(Q_down, Q, rf_h, Q_state, rf_decim);
        
        //Demodulation
        fmDemodulation(I_down, Q_down, prev_I, prev_Q,IF);
        
        
		
        // --------------- MONO PATH ---------------------------
		
		//impulseResponseLPF(audio_IF,audio_Fc,rf_taps,audio_h);
		
		//auto start3 = std::chrono::high_resolution_clock::now();
        //Demodulation
        impulseResponseLPF(audio_IF,audio_Fc,rf_taps, audio_h);
        decim4(IF_down, IF, audio_h, IF_state, audio_decim, audio_up);
        //auto stop3 = std::chrono::high_resolution_clock::now();
        
        //auto duration3 =  std::chrono::duration<double, std::milli>(stop3 - start3);

		// Output time taken
		//std::cerr << "Time taken by upsampling function: "
              //<< duration3.count() << " microseconds" << std::endl;
		
		/*
		// --------------- STEREO PATH ----------------------
		// Extraction

		
		decim2(extraction_filtered, IF, extraction_h, extraction_state, 1);
		//std::cerr << "size of extraction filtered: " << extraction_filtered.size() << std::endl;

		
		decim2(carrier_filtered, IF, carrier_h, carrier_state, 1);
		//std::cerr << "size of carrier filtered: " << carrier_filtered.size() << std::endl;
		
		// PLL and nco
		fmPLL(carrier_filtered, PLL_frequency, audio_IF, nco_out, integrator, phase_estimate, feedbackI, feedbackQ, trig_offset, nco_scale, phase_adjust, norm_bandwith);
		
		// Mixing
        for (int i = 0; i < extraction_filtered.size(); i++) {
                mixed_stereo[i] = extraction_filtered[i]*nco_out[i]*2;
        }
        //std::cerr << "size of mixed_stereo: " << mixed_stereo.size() << std::endl;
            
			
		
		decim4(IF_down, IF, audio_h, IF_state, audio_decim, audio_up);
		//std::cerr << "size of IF DOWN: " << IF_down.size() << std::endl;
		//std::cerr << "size of IF DOWN: " << IF_down.size() << std::endl;
		delayBlock(IF_down, IF_delayed, delay_state);
		//std::cerr << "size of extraction: " << extraction_delayed.size() << std::endl;
			

		decim4(stereo_conv, mixed_stereo, stereo_h, stereo_state, audio_decim, audio_up);	
		//std::cerr << "size of stereo conv: " << stereo_conv.size() << std::endl;
			

		// STEREO COMBINATION
		for (int j = 0; j < stereo_conv.size(); j++){
			stereo_left_right[2*j] = (IF_delayed[j] + stereo_conv[j])/2;
			stereo_left_right[2*j+1] = (IF_delayed[j] - stereo_conv[j])/2;
		}
		//std::cerr << "size of stereo LR: " << stereo_conv.size() << std::endl;
		
		
		/*
		std::cerr << "Size of BPF coeff: " << carrier_h.size() << std::endl;
		std::cerr << "carrier_h: ";
		for (unsigned int i = 0; i <  carrier_h.size(); i++) {
			std::cerr << carrier_h[i] << ", ";
		}
		std::cerr << std::endl;
		*/
		
		//convolveFIR_state(carrier_filtered, IF, carrier_h, carrier_state);
		
		//impulseResponseBPF(audio_IF, extraction_fb, extraction_fe, n_taps, extraction_h);
		//convolveFIR_state(extraction_filtered, IF, extraction_h, extraction_state);
		
		
	    /*
		//FOR TESTING PURPOSES
		std::cerr << "Size of block_data: " << block_data.size() << std::endl;
		std::cerr << "Size of I: " << I.size() << std::endl;
		std::cerr << "Size of Q: " << Q.size() << std::endl;
		std::cerr << "Size of rf_h: " << rf_h.size() << std::endl;
		std::cerr << "Size of I_filtered: " << I_filtered.size() << std::endl;
		std::cerr << "Size of Q_filtered: " << Q_filtered.size() << std::endl;
		std::cerr << "Size of I_down: " << I_down.size() << std::endl;
		std::cerr << "Size of Q_down: " << Q_down.size() << std::endl;
		std::cerr << "Size of IF: " << IF.size() << std::endl;
		std::cerr << "Size of audio_h: " << audio_h.size() << std::endl;
		std::cerr << "Size of IF_filtered: " << IF_filtered.size() << std::endl;
		std::cerr << "Size of IF_down: " << IF_down.size() << std::endl;
		

        std::cerr << "IF: ";
		for (unsigned int i = 0; i <  IF.size(); i++) {
			std::cerr << IF[i] << " ";
		}
		std::cerr << std::endl; 
		*/


        //-------------------------------------------------------------------------------WRITING TO APLAY------------------------------------------------------------------------------------------------------------------
        
		std::vector<short int> audio_data(IF_down.size());
		for (unsigned int k=0; k<IF_down.size(); k++) {
			if (std::isnan(IF_down[k])) audio_data[k] = 0;
			// prepare a block of audio data to be redirected to stdout at once
			else audio_data[k] = static_cast<short int>(IF_down[k] * 16384);
		}

		// a block write is approx an order of magnitude faster than writing each sample
		// (assuming input streams are processed in blocks of hundreds of kilobytes)
		fwrite(&audio_data[0], sizeof(short int), audio_data.size(), stdout);  
		

        //------------------------------------------------------------------------------- done? --------------------------------------------------------------
  
  }

	return 0;
}
