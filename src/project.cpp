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


int main(int argc, char* argv[])
{
	// ------------------------------------------------------------MODE SELECTION---------------------------------------------------------------------
		// using namespace std::chrono;
		//default mode 0 
		int mode = 0;
		std::string monoStereoMode;
		
		if (argc<2){
			std::cerr << "Operating in default mode 0 with Mono" << std::endl;
			
		} else if (argc==2){
			mode=atoi(argv[1]);
			
			if (mode>3){
				std::cerr << "Wrong mode " << mode << std::endl;
				exit(1);
			}
        std::cerr << "Operating in mode " << mode << " with Mono output" << std::endl;
		} else if (argc ==3) {
			mode = atoi(argv[1]);
			monoStereoMode = (argv[2]);
			
			if (mode > 3 || (monoStereoMode != "m" && monoStereoMode != "s")) {
				std::cerr << "Invalid arguments. Mode should be between 0-3 and Mono/Stereo should be 'm' or 's'." << std::endl;
				exit(1);
			}
			
			std::cerr << "Operating in mode " << mode << " with " << (monoStereoMode == "m" ? "Mono" : "Stereo") << " output" << std::endl;
		} else {
			std::cerr << "Usage: " << argv[0] << " <mode> <Mono/Stereo ('m'/'s')>" << std::endl;
			exit(1);
		}
       
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
		
		float n_taps = 101;
		
		// -------------------- MONO ------------------------
		
		float rf_Fc = 100e3;
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
        std::vector<float> IF_state(n_taps - 1, 0.0);
        
		std::vector<float> I_state(n_taps - 1, 0.0);
        std::vector<float> Q_state(n_taps - 1, 0.0); 
        
        std::vector<float> I_filtered;
        std::vector<float> Q_filtered;
        
        std::vector<float> mono_output;
       
        // -------------------- STEREO -----------------------    

        std::vector<float> extraction_h;
        std::vector<float> extraction_filtered;
        std::vector<float> carrier_h;
        std::vector<float> carrier_filtered;
        std::vector<float> carrier_state(n_taps - 1, 0.0);
        std::vector<float> extraction_state(n_taps - 1, 0.0);
        
        float extraction_fb = 22e3;
        float extraction_fe = 54e3;
        float carrier_fb = 18.5e3;
        float carrier_fe = 19.5e3;
        
        std::vector<float> mixed_stereo(BLOCK_SIZE/(2*rf_decim), 0.0);
        
        std::vector<float> stereo_state(n_taps, 0.0);
        std::vector<float> stereo_conv;
		std::vector<float> stereo_h;
		std::vector<float> stereo_left_right;
		
		float PLL_frequency = 19e3;
		float integrator = 0.0;
		float phase_estimate = 0.0;
		float feedbackI = 1.0;
		float feedbackQ = 0.0; 
		float trig_offset = 0.0; 
		float nco_scale = 2.0;
		float phase_adjust = 0.0; 
		float norm_bandwith = 0.01;
		std::vector<float> nco_out;
		std::vector<float> pll_state = {integrator, phase_estimate, feedbackI, feedbackQ, trig_offset, 1.0};
	
        std::vector<float> IF_delayed;
        std::vector<float> delay_state((n_taps - 1)/2, 0.0);
        
        // -------------------- RDS -----------------------
        
      
        float rds_extraction_fb = 54e3;
        float rds_extraction_fe = 60e3;
        std::vector<float> rds_extraction_h;
        float rds_carrier_fb = 113.5e3;
        float rds_carrier_fe = 114.5e3;
        std::vector<float> rds_carrier_h;
        
        std::vector<float> rds_extraction_state(n_taps - 1, 0.0);
        std::vector<float> rds_extraction_filtered;
        std::vector<float> rds_carrier_state(n_taps - 1, 0.0);
        std::vector<float> rds_carrier_filtered;
        std::vector<float> rds_squared;
        
        std::vector<float> rds_delayed;
        std::vector<float> rds_delay_state((n_taps - 1)/2, 0.0);
        
        float rds_PLL_frequency = 114e3;
		float rds_integrator = 0.0;
		float rds_phase_estimate = 0.0;
		float rds_feedbackI = 1.0;
		float rds_feedbackQ = 0.0; 
		float rds_trig_offset = 0.0; 
		float rds_nco_scale = 2.0;
		float rds_phase_adjust = 0.0; 
		float rds_norm_bandwith = 0.01;
		std::vector<float> rds_nco_out;
		std::vector<float> rds_pll_state = {rds_integrator, rds_phase_estimate, rds_feedbackI, rds_feedbackQ, rds_trig_offset, 1.0};
		
		std::vector<float> mixed_audio(BLOCK_SIZE/(2*rf_decim), 0.0);
		
		
        
        // -------------------- FILTER COEFFICIENTS -----------------------  
        
        impulseResponseLPF(rf_Fs,rf_Fc,n_taps,rf_h);
        impulseResponseBPF(audio_IF, extraction_fb, extraction_fe, n_taps, extraction_h);
        impulseResponseBPF(audio_IF, carrier_fb, carrier_fe, n_taps, carrier_h);
        impulseResponseLPF(audio_IF,audio_Fc,n_taps, audio_h);
        impulseResponseLPF(audio_IF, audio_Fc, n_taps, stereo_h);
           
        
        // ------------------ BLOCK READING ----------------------
        
       for(unsigned int block_id = 0; ; block_id++){
		
		//Read raw data 
		readStdinBlockData(BLOCK_SIZE, block_id, block_data);
		
		if((std::cin.rdstate()) != 0){
            std::cerr << "End of input stream reached!" << std::endl;
            exit(1);
        }
        
        std::cerr << "Processing block " << block_id <<std::endl;  //output  message indicating that a block has been successfully read
        
        
		// ------------------- RF FRONT END ---------------------
		// Splitting block data -> I, Q
		
        int iq_index = 0;
        for (unsigned int i=0; i < block_data.size(); i+=2) {
            I[iq_index] = block_data[i];
            Q[iq_index] = block_data[i + 1];
            iq_index++;
        }
        
        //Decimation
        decim2(I_down, I, rf_h, I_state, rf_decim);       
        decim2(Q_down, Q, rf_h, Q_state, rf_decim);
        
        //Demodulation
        fmDemodulation(I_down, Q_down, prev_I, prev_Q,IF);
        
		
        // ------------------ MONO PATH --------------------------
		
		//impulseResponseLPF(audio_IF,audio_Fc,rf_taps,audio_h);
		
		//auto start3 = std::chrono::high_resolution_clock::now();
        //Demodulation
        //impulseResponseLPF(audio_IF,audio_Fc,rf_taps, audio_h);
        //decim4(IF_down, IF, audio_h, IF_state, audio_decim, audio_up);
        //auto stop3 = std::chrono::high_resolution_clock::now();
        
        //auto duration3 =  std::chrono::duration<double, std::milli>(stop3 - start3);

		// Output time taken
		//std::cerr << "Time taken by upsampling function: "
              //<< duration3.count() << " microseconds" << std::endl;
		
		
		// --------------- STEREO PATH ----------------------
		
		// carrier path
		decim2(carrier_filtered, IF, carrier_h, carrier_state, 1);
		fmPLL(carrier_filtered, PLL_frequency, audio_IF, nco_scale, phase_adjust, norm_bandwith, pll_state, nco_out);
		
		// stereo channel extraction
		decim2(extraction_filtered, IF, extraction_h, extraction_state, 1);
		
		// stereo processing
		stereo_mixer(mixed_stereo, extraction_filtered, nco_out);
		decim4(stereo_conv, mixed_stereo, stereo_h, stereo_state, audio_decim, audio_up);
		
		// delay matching
		all_pass_filter(IF_delayed, IF, delay_state, n_taps);
		decim4(IF_down, IF_delayed, audio_h, IF_state, audio_decim, audio_up);
		
		// combiner for left/right channel audio
		stereo_combiner(stereo_left_right, IF_down, stereo_conv);
		
		// for testing mono
		//decim4(IF_down, IF, audio_h, IF_state, audio_decim, audio_up);
		
		// ------------------- RDS PATH ----------------------
        /*
        // filter coefficients
        impulseResponseBPF(audio_IF, rds_extraction_fb, rds_extraction_fe, n_taps, rds_extraction_h);
        impulseResponseBPF(audio_IF, rds_carrier_fb, rds_carrier_fe, n_taps, rds_carrier_h);
        
        // rds channel
        decim2(rds_extraction_filtered, IF, rds_extraction_h, rds_extraction_state, 1);
        all_pass_filter(rds_delayed, rds_extraction_filtered, rds_delay_state, n_taps);
        
        // rds carrier
        squaring_nonlinearity(rds_extraction_filtered, rds_squared)
        decim2(rds_carrier_filtered, rds_squared, rds_carrier_h, rds_carrier_state, 1);
		fmPLL(rds_carrier_filtered, rds_PLL_frequency, audio_IF, rds_nco_scale, rds_phase_adjust, rds_norm_bandwith, rds_pll_state, rds_nco_out);
		
		// mixing
        stereo_mixer(mixed_audio, rds_delayed, rds_nco_out);
        */
		
		// ------------------- TESTING ----------------------
		
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
        if (monoStereoMode == "m"){
		std::vector<short int> audio_data(IF_down.size());
		for (unsigned int k=0; k<IF_down.size(); k++) {
			if (std::isnan(IF_down[k])) audio_data[k] = 0;
			// prepare a block of audio data to be redirected to stdout at once
			else audio_data[k] = static_cast<short int>(IF_down[k] * 16384);
			}
			fwrite(&audio_data[0], sizeof(short int), audio_data.size(), stdout); 
		}
		
        else{
		std::vector<short int> audio_data(stereo_left_right.size());
		for (unsigned int k=0; k<stereo_left_right.size(); k++) {
			if (std::isnan(stereo_left_right[k])) audio_data[k] = 0;
			// prepare a block of audio data to be redirected to stdout at once
			else audio_data[k] = static_cast<short int>(stereo_left_right[k] * 16384);
			}
			fwrite(&audio_data[0], sizeof(short int), audio_data.size(), stdout); 
		// a block write is approx an order of magnitude faster than writing each sample
		// (assuming input streams are processed in blocks of hundreds of kilobytes)
		}
		 
        //------------------------------------------------------------------------------- done? --------------------------------------------------------------
  
  }

	return 0;
}
