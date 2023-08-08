%% For use in Simulink Model

%% Notes
%{
    - If the start frequency is higher than 77 GHz, there is a little bit
    of extra work that has to be done. When creating this simulation, I
    removed this functionality to simplify things. I did implement the
    functionality in previous versions though which can be found in the
    archived functions folder. See the configure_FMCW_parameters function
    for this information

    - There are a lot of other default radar configuration scenarios that I
    have computed, but in this simulation I'm only using the "realistic"
    settings for now. The other configurations can be found in the archived
    Simulator code

    - Chirp construction is slightly different. See the notes in the Radar
    class for more details

    - Removed all of the plotting functionality in this revision as it
    won't be used in these simulations
%}

classdef Simulator_revA < handle
    %SIMULATOR class used to support simulation in the simulink model

    properties (Access = public)
        Attacker                % an Attacker() class object
        Victim                  % a Radar() class object
        SimulatedTarget         % a Target() class object 

        %variables used in computing the FMCW simulation - MOVE TO PRIVATE
        channel_target          %phased.FreeSpace object for a target
        channel_attacker        %phased.FreeSpace object for an attacker

        %the following properties are used as support parameters when
        %simulating the interaction between various parts of the FMCW
        %simulation
        sensing_subsystem_support
    end

    %the following properties are used as support parameters when
    %simulating the interaction between various parts of the FMCW
    %simulation
    properties (Access = private)
        
    end
    
    methods (Access = public)
        function obj = Simulator_revA()
            %Simulator Construct an instance of this class
            %   Detailed explanation goes here
            obj.Attacker = Attacker_revA();
            obj.Victim = Radar_revA();
        end

%% [1] Functions to configure the simulation
        
        function configure_FMCW_Radar_parameters(obj)
            %{
                Purpose: initializes all of the needed waveform parameters
                    to simulate the actual waveforms for the attacker and
                    victim
                Note: The attacker and Victim radar components will have
                the same setting. In practice this might not be the case,
                but for the sake of simplifying the simulation, they will
                have the same settings
            %}

            %obj.Victim.downsample_factor = ceil((2 * obj.Victim.Chirp_Tx_Bandwidth_MHz * 1e6) / (obj.Victim.ADC_SampleRate_MSps * 1e6));
                %NOTE: if needed, can use nextpow2() function when
                %implementing on the fpga

            obj.Victim.downsample_factor = ceil((obj.Victim.Chirp_Tx_Bandwidth_MHz * 1e6) / (obj.Victim.ADC_SampleRate_MSps * 1e6));
            obj.Victim.Radar_Signal_Processor.decimation_factor = obj.Victim.downsample_factor;
            obj.Attacker.Subsystem_tracking.downsample_factor = obj.Victim.downsample_factor;
            obj.Attacker.Subsystem_tracking.Radar_Signal_Processor.decimation_factor = obj.Victim.downsample_factor;
            
            %removed the doubling of the chirp Tx Bandwidth since we are
            %assuming complex sampling

            obj.Victim.FMCW_sampling_rate_Hz = obj.Victim.ADC_SampleRate_MSps * 1e6 * obj.Victim.downsample_factor;
            obj.Attacker.Subsystem_tracking.FMCW_sampling_rate_Hz = obj.Victim.FMCW_sampling_rate_Hz;

            %set the FMCW sampling period as it is good to have for
            %reference
            obj.Victim.FMCW_sampling_period_s = 1/obj.Victim.FMCW_sampling_rate_Hz;
            obj.Attacker.Subsystem_tracking.FMCW_sampling_period_s = obj.Victim.FMCW_sampling_period_s;
            

            %set sweep time
            obj.Victim.sweep_time = obj.Victim.RampEndTime_us * 1e-6;
            obj.Attacker.Subsystem_tracking.sweep_time = obj.Attacker.Subsystem_tracking.RampEndTime_us * 1e-6;
            
            
            %configure the FMCW waveforms for the victim and attacker
            %configure waveforms
            obj.Victim.configure_waveform_and_chirp();
            obj.Attacker.Subsystem_tracking.configure_waveform_and_chirp();
            
            obj.channel_target = phased.FreeSpace( ...
                "PropagationSpeed",physconst('LightSpeed'), ...
                "OperatingFrequency", obj.Victim.StartFrequency_GHz * 1e9, ...
                "SampleRate", obj.Victim.FMCW_sampling_rate_Hz, ...
                "TwoWayPropagation", true);

            obj.channel_attacker = phased.FreeSpace( ...
                "PropagationSpeed",physconst('LightSpeed'), ...
                "OperatingFrequency", obj.Victim.StartFrequency_GHz * 1e9, ...
                "SampleRate", obj.Victim.FMCW_sampling_rate_Hz, ...
                "TwoWayPropagation", false);
            
            
            %configure transmitters, receivers, lowpass filters, and CFAR
            %detectors

            obj.Attacker.Subsystem_tracking.configure_transmitter_and_receiver();
            obj.Attacker.Subsystem_tracking.configure_radar_signal_processor();

            obj.Victim.configure_transmitter_and_receiver();
            obj.Victim.configure_radar_signal_processor();
            %
        end

        function load_params_from_JSON(obj,file_path)
            fileID = fopen(file_path,"r");
            json_text = fileread(fileID);
            json = jsondecode(json_text);
        end
        
        function load_realistic_victim_params(obj)
            %setup the victim's chirp parameters
            obj.Victim.StartFrequency_GHz         = 77.0;
            obj.Victim.FrequencySlope_MHz_us      = 10.76;
            obj.Victim.TxStartTime_us             = 0;
            obj.Victim.ADC_Samples                = 256;
            obj.Victim.ADC_SampleRate_MSps        = 7.17;
            obj.Victim.ChirpCycleTime_us          = 50;             

            %setup the victim's frame parameters
            obj.Victim.NumChirps                  = 128;
            obj.Victim.FramePeriodicity_ms        = 33.33;
            
            %define plot color default values
            obj.Victim.plotResolution_us = 0.01;
            obj.Victim.tx_period_plot_color = 'blue';
            obj.Victim.tx_sampling_period_plot_color = 'cyan';
            obj.Victim.radar_name = 'Victim';

            %set timing offset to zero as this is the victim
            obj.Victim.timing_offset_us = 0;

            %compute all remaining "calculated" values
            obj.Victim.compute_calculated_vals();
        end

        function load_realisitc_attacker_params(obj)
            %setup the victim's chirp parameters
            obj.Attacker.Subsystem_tracking.StartFrequency_GHz         = 77.0;
            obj.Attacker.Subsystem_tracking.FrequencySlope_MHz_us      = 10.76;
            obj.Attacker.Subsystem_tracking.TxStartTime_us             = 0;
            obj.Attacker.Subsystem_tracking.ADC_Samples                = 256;
            obj.Attacker.Subsystem_tracking.ADC_SampleRate_MSps        = 7.17;
            obj.Attacker.Subsystem_tracking.ChirpCycleTime_us          = 50;
            
            %setup the victim's frame parameters
            obj.Attacker.Subsystem_tracking.NumChirps                  = 128;
            obj.Attacker.Subsystem_tracking.FramePeriodicity_ms        = 33.4965;
            
            %define plot color default values
            obj.Attacker.Subsystem_tracking.plotResolution_us = 0.01;
            obj.Attacker.Subsystem_tracking.tx_period_plot_color = 'red';
            obj.Attacker.Subsystem_tracking.tx_sampling_period_plot_color = 'magenta';
            obj.Attacker.Subsystem_tracking.radar_name = 'Attacker';

            %define set the default offset to be 0us
            obj.Attacker.Subsystem_tracking.timing_offset_us = 0;

            %compute all remaining "calculated" values
            obj.Attacker.Subsystem_tracking.compute_calculated_vals();
        end

        function load_B210_victim_params(obj)
            %setup the victim's chirp parameters
            obj.Victim.StartFrequency_GHz         = 77.0;
            obj.Victim.FrequencySlope_MHz_us      = 1.2;
            obj.Victim.TxStartTime_us             = 0;
            obj.Victim.ADC_Samples                = 64;
            obj.Victim.ADC_SampleRate_MSps        = 1.54;
            obj.Victim.ChirpCycleTime_us          = 49.98;             

            %setup the victim's frame parameters
            obj.Victim.NumChirps                  = 64;
            obj.Victim.FramePeriodicity_ms        = 33.33;
            
            %define plot color default values
            obj.Victim.plotResolution_us = 0.01;
            obj.Victim.tx_period_plot_color = 'blue';
            obj.Victim.tx_sampling_period_plot_color = 'cyan';
            obj.Victim.radar_name = 'Victim';

            %set timing offset to zero as this is the victim
            obj.Victim.timing_offset_us = 0;

            %compute all remaining "calculated" values
            obj.Victim.compute_calculated_vals();
        end

        function load_B210_victim_params_highvres(obj)
            %setup the victim's chirp parameters
            obj.Victim.StartFrequency_GHz         = 5.8;
            obj.Victim.FrequencySlope_MHz_us      = 1;
            obj.Victim.TxStartTime_us             = 0;
            obj.Victim.ADC_Samples                = 64;
            obj.Victim.ADC_SampleRate_MSps        = 3.2;
            obj.Victim.ChirpCycleTime_us          = 600;             

            %setup the victim's frame parameters
            obj.Victim.NumChirps                  = 128;
            obj.Victim.FramePeriodicity_ms        = 200;
            
            %define plot color default values
            obj.Victim.plotResolution_us = 0.01;
            obj.Victim.tx_period_plot_color = 'blue';
            obj.Victim.tx_sampling_period_plot_color = 'cyan';
            obj.Victim.radar_name = 'Victim';

            %set timing offset to zero as this is the victim
            obj.Victim.timing_offset_us = 0;

            %compute all remaining "calculated" values
            obj.Victim.compute_calculated_vals();
        end
        
        function load_B210_victim_params_lowBW(obj)
            %setup the victim's chirp parameters
            obj.Victim.StartFrequency_GHz         = 2;
            obj.Victim.FrequencySlope_MHz_us      = 1;
            obj.Victim.TxStartTime_us             = 0;
            obj.Victim.ADC_Samples                = 64;
            obj.Victim.ADC_SampleRate_MSps        = 3.2;
            obj.Victim.ChirpCycleTime_us          = 45;             

            %setup the victim's frame parameters
            obj.Victim.NumChirps                  = 256;
            obj.Victim.FramePeriodicity_ms        = 33.33;
            
            %define plot color default values
            obj.Victim.plotResolution_us = 0.01;
            obj.Victim.tx_period_plot_color = 'blue';
            obj.Victim.tx_sampling_period_plot_color = 'cyan';
            obj.Victim.radar_name = 'Victim';

            %set timing offset to zero as this is the victim
            obj.Victim.timing_offset_us = 0;

            %compute all remaining "calculated" values
            obj.Victim.compute_calculated_vals();
        end

        function load_B210_victim_params_highBW(obj)
            %setup the victim's chirp parameters
            obj.Victim.StartFrequency_GHz         = 5.8;
            obj.Victim.FrequencySlope_MHz_us      = 1;
            obj.Victim.TxStartTime_us             = 0;
            obj.Victim.ADC_Samples                = 64;
            obj.Victim.ADC_SampleRate_MSps        = 1.83;
            obj.Victim.ChirpCycleTime_us          = 50;             

            %setup the victim's frame parameters
            obj.Victim.NumChirps                  = 256;
            obj.Victim.FramePeriodicity_ms        = 33.33;
            
            %define plot color default values
            obj.Victim.plotResolution_us = 0.01;
            obj.Victim.tx_period_plot_color = 'blue';
            obj.Victim.tx_sampling_period_plot_color = 'cyan';
            obj.Victim.radar_name = 'Victim';

            %set timing offset to zero as this is the victim
            obj.Victim.timing_offset_us = 0;

            %compute all remaining "calculated" values
            obj.Victim.compute_calculated_vals();
        end

        function load_B210_victim_params_100MHzBW(obj)
            %setup the victim's chirp parameters
            obj.Victim.StartFrequency_GHz         = 1.5;
            obj.Victim.FrequencySlope_MHz_us      = 2.13;
            obj.Victim.TxStartTime_us             = 0;
            obj.Victim.ADC_Samples                = 64;
            obj.Victim.ADC_SampleRate_MSps        = 1.6;
            obj.Victim.ChirpCycleTime_us          = 50;             

            %setup the victim's frame parameters
            obj.Victim.NumChirps                  = 256;
            obj.Victim.FramePeriodicity_ms        = 33.33;
            
            %define plot color default values
            obj.Victim.plotResolution_us = 0.01;
            obj.Victim.tx_period_plot_color = 'blue';
            obj.Victim.tx_sampling_period_plot_color = 'cyan';
            obj.Victim.radar_name = 'Victim';

            %set timing offset to zero as this is the victim
            obj.Victim.timing_offset_us = 0;

            %compute all remaining "calculated" values
            obj.Victim.compute_calculated_vals();
        end
        
        function load_B210_attacker_params(obj)
            %setup the victim's chirp parameters
            obj.Attacker.Subsystem_tracking.StartFrequency_GHz                  = 77.0;
            obj.Attacker.Subsystem_tracking.FrequencySlope_MHz_us               = 1.2;
            obj.Attacker.Subsystem_tracking.TxStartTime_us                      = 0;
            obj.Attacker.Subsystem_tracking.ADC_Samples                         = 64;
            obj.Attacker.Subsystem_tracking.ADC_SampleRate_MSps                 = 1.54;
            obj.Attacker.Subsystem_tracking.ChirpCycleTime_us                   = 49.98;             

            %setup the victim's frame parameters
            obj.Attacker.Subsystem_tracking.NumChirps                  = 64;
            obj.Attacker.Subsystem_tracking.FramePeriodicity_ms        = 33.33;
            
            %define plot color default values
            obj.Attacker.Subsystem_tracking.plotResolution_us = 0.01;
            obj.Attacker.Subsystem_tracking.tx_period_plot_color = 'red';
            obj.Attacker.Subsystem_tracking.tx_sampling_period_plot_color = 'magenta';
            obj.Attacker.Subsystem_tracking.radar_name = 'Attacker';

            %set timing offset to zero as this is the victim
            obj.Attacker.Subsystem_tracking.timing_offset_us = 0;

            %compute all remaining "calculated" values
            obj.Attacker.Subsystem_tracking.compute_calculated_vals();
        end
        
        function load_B210_attacker_params_lowBW(obj)
            %setup the attacker's chirp parameters
            obj.Attacker.Subsystem_tracking.StartFrequency_GHz         = 2;
            obj.Attacker.Subsystem_tracking.FrequencySlope_MHz_us      = 1;
            obj.Attacker.Subsystem_tracking.TxStartTime_us             = 0;
            obj.Attacker.Subsystem_tracking.ADC_Samples                = 64;
            obj.Attacker.Subsystem_tracking.ADC_SampleRate_MSps        = 3.2;
            obj.Attacker.Subsystem_tracking.ChirpCycleTime_us          = 45;             

            %setup the attacker's frame parameters
            obj.Attacker.Subsystem_tracking.NumChirps                  = 256;
            obj.Attacker.Subsystem_tracking.FramePeriodicity_ms        = 33.33;
            
            %define plot color default values
            obj.Attacker.Subsystem_tracking.plotResolution_us = 0.01;
            obj.Attacker.Subsystem_tracking.tx_period_plot_color = 'blue';
            obj.Attacker.Subsystem_tracking.tx_sampling_period_plot_color = 'cyan';
            obj.Attacker.Subsystem_tracking.radar_name = 'Victim';

            %set timing offset to zero as this is the victim
            obj.Attacker.Subsystem_tracking.timing_offset_us = 0;

            %compute all remaining "calculated" values
            obj.Attacker.Subsystem_tracking.compute_calculated_vals();
        end

        function load_B210_attacker_params_highBW(obj)
            %setup the attacker's chirp parameters
            obj.Attacker.Subsystem_tracking.StartFrequency_GHz         = 5.8;
            obj.Attacker.Subsystem_tracking.FrequencySlope_MHz_us      = 1;
            obj.Attacker.Subsystem_tracking.TxStartTime_us             = 0;
            obj.Attacker.Subsystem_tracking.ADC_Samples                = 64;
            obj.Attacker.Subsystem_tracking.ADC_SampleRate_MSps        = 1.83;
            obj.Attacker.Subsystem_tracking.ChirpCycleTime_us          = 50;             

            %setup the attacker's frame parameters
            obj.Attacker.Subsystem_tracking.NumChirps                  = 256;
            obj.Attacker.Subsystem_tracking.FramePeriodicity_ms        = 33.33;
            
            %define plot color default values
            obj.Attacker.Subsystem_tracking.plotResolution_us = 0.01;
            obj.Attacker.Subsystem_tracking.tx_period_plot_color = 'blue';
            obj.Attacker.Subsystem_tracking.tx_sampling_period_plot_color = 'cyan';
            obj.Attacker.Subsystem_tracking.radar_name = 'Victim';

            %set timing offset to zero as this is the victim
            obj.Attacker.Subsystem_tracking.timing_offset_us = 0;

            %compute all remaining "calculated" values
            obj.Attacker.Subsystem_tracking.compute_calculated_vals();
        end

        function load_B210_attacker_params_100MHzBW(obj)
            %setup the attacker's chirp parameters
            obj.Attacker.Subsystem_tracking.StartFrequency_GHz         = 1.5;
            obj.Attacker.Subsystem_tracking.FrequencySlope_MHz_us      = 2.13;
            obj.Attacker.Subsystem_tracking.TxStartTime_us             = 0;
            obj.Attacker.Subsystem_tracking.ADC_Samples                = 64;
            obj.Attacker.Subsystem_tracking.ADC_SampleRate_MSps        = 1.6;
            obj.Attacker.Subsystem_tracking.ChirpCycleTime_us          = 50;             

            %setup the attacker's frame parameters
            obj.Attacker.Subsystem_tracking.NumChirps                  = 256;
            obj.Attacker.Subsystem_tracking.FramePeriodicity_ms        = 33.33;
            
            %define plot color default values
            obj.Attacker.Subsystem_tracking.plotResolution_us = 0.01;
            obj.Attacker.Subsystem_tracking.tx_period_plot_color = 'blue';
            obj.Attacker.Subsystem_tracking.tx_sampling_period_plot_color = 'cyan';
            obj.Attacker.Subsystem_tracking.radar_name = 'Victim';

            %set timing offset to zero as this is the victim
            obj.Attacker.Subsystem_tracking.timing_offset_us = 0;

            %compute all remaining "calculated" values
            obj.Attacker.Subsystem_tracking.compute_calculated_vals();
        end

        function load_B210_attacker_params_highvres(obj)
            %setup the attacker's chirp parameters
            obj.Attacker.Subsystem_tracking.StartFrequency_GHz         = 5.8;
            obj.Attacker.Subsystem_tracking.FrequencySlope_MHz_us      = 1;
            obj.Attacker.Subsystem_tracking.TxStartTime_us             = 0;
            obj.Attacker.Subsystem_tracking.ADC_Samples                = 64;
            obj.Attacker.Subsystem_tracking.ADC_SampleRate_MSps        = 3.2;
            obj.Attacker.Subsystem_tracking.ChirpCycleTime_us          = 600;             

            %setup the attacker's frame parameters
            obj.Attacker.Subsystem_tracking.NumChirps                  = 128;
            obj.Attacker.Subsystem_tracking.FramePeriodicity_ms        = 200;
            
            %define plot color default values
            obj.Attacker.Subsystem_tracking.plotResolution_us = 0.01;
            obj.Attacker.Subsystem_tracking.tx_period_plot_color = 'blue';
            obj.Attacker.Subsystem_tracking.tx_sampling_period_plot_color = 'cyan';
            obj.Attacker.Subsystem_tracking.radar_name = 'Victim';

            %set timing offset to zero as this is the victim
            obj.Attacker.Subsystem_tracking.timing_offset_us = 0;

            %compute all remaining "calculated" values
            obj.Attacker.Subsystem_tracking.compute_calculated_vals();
        end

        function load_target_realistic(obj)
            %{
                Purpose: load a realistic target
            %}

            %configure the simulated target
            position_m = [55;0;0];
            velocity_meters_per_s = [-5;0;0];
            rcs_sq_meters = db2pow(min(10*log10(norm(position_m))+5,20));
            operating_frequency_Hz = obj.Victim.StartFrequency_GHz * 1e9;

            obj.SimulatedTarget = Target_revA(position_m,velocity_meters_per_s,rcs_sq_meters,operating_frequency_Hz);
        end
        
        function load_realistic_attacker_and_victim_position_and_velocity(obj)
            %{
                Purpose: configures a default scenario for the attacker and
                victim positions and velocities
            %}
            obj.Victim.position_m = [0;0;0];
            obj.Victim.velocity_m_per_s = [7;0;0];
            obj.Victim.platform = phased.Platform( ...
                'InitialPosition',obj.Victim.position_m, ...
                'Velocity',obj.Victim.velocity_m_per_s);

            obj.Attacker.position_m = [75;0;0];
            obj.Attacker.velocity_m_per_s = [0;0;0];
            obj.Attacker.configure_platform();
        end
    
        function initialize_sensing_subsystem_support(obj)
            %{
                Purpose: the purpose of this function is to initialize
                private parameters that will be used in assembling the
                signal that the sensing subsystem will receive from victim
            %}
            %parameters to keep track of the index from the transmitted
            %chirp signal from the victim
            obj.sensing_subsystem_support.chirp_sample_index = 1;
            obj.sensing_subsystem_support.max_chirp_sample_index = size(obj.Victim.chirp,1);
            obj.sensing_subsystem_support.chirp_samples_left = obj.sensing_subsystem_support.max_chirp_sample_index;
            
            %parameters to keep track of the index from the transmitted
            %signal frame from the victim
            obj.sensing_subsystem_support.max_frame_sample_index = obj.Victim.num_samples_per_frame;
            obj.sensing_subsystem_support.frame_samples_left = obj.sensing_subsystem_support.max_frame_sample_index;
             
            %parameters to keep track of the index of the received signal
            %received by the sensing subsystem
            obj.sensing_subsystem_support.received_sample_index = 1;
            obj.sensing_subsystem_support.max_received_sample_index = obj.Attacker.Subsystem_spectrum_sensing.spectogram_params.num_ADC_samples_per_spectogram;
            obj.sensing_subsystem_support.received_samples_left = obj.sensing_subsystem_support.max_received_sample_index;
            
            %parameters to keep track of the index of how many chirps and
            %frames the victim has transmitted
            obj.sensing_subsystem_support.chirp_count = 0;
            obj.sensing_subsystem_support.frame_count = 0;
        end

%% [2] Functions for running the FMCW Simulation on Matlab

        function [victim_pos, victim_vel,attacker_pos, attacker_vel, tgt_pos,tgt_vel] = FMCW_determine_positions_and_velocities(obj,victim_frame,victim_chirp)
        %{
            Purpose: determine the position and velocity for the
                attacker,defender, and target at the end of a given chirp
                in a given frame
            Inputs:
                victim_chirp: the desired chirp
                victim_frame: the desired frame
            Outputs:
                [victim_pos, victim_vel]: the position and velocity of
                    the victim at the start of the desired chirp
                [attacker_pos, attacker_vel]: the position and velocity of
                    the attacker at the start of the desired chirp
                [target_pos, target_vel]: the position and velocity of
                    the attacker at the start of the desired chirp
        %}
            %calculate the time that the desired chirp will end at
            frame_start_time_s = (obj.Victim.FramePeriodicity_ms * 1e-3) * (victim_frame - 1);
            chirp_end_time = frame_start_time_s + (obj.Victim.ChirpCycleTime_us * 1e-6) *...
                (victim_chirp);
            
            %reset the positions of each platform to be safe
            reset(obj.Victim.platform);
            reset(obj.Attacker.platform);
            reset(obj.SimulatedTarget.platform);

            if chirp_end_time <= 0
                time_increment = obj.Victim.ChirpCycleTime_us * 1e-6;
                [victim_pos, victim_vel] = obj.Victim.platform(time_increment);
                [tgt_pos,tgt_vel] = obj.SimulatedTarget.platform(time_increment);
                [attacker_pos, attacker_vel] = obj.Attacker.platform(time_increment);
            else
                %take a step for each platform so that they update to the
                %values at the start of the frame
                obj.Victim.platform(chirp_end_time);
                obj.SimulatedTarget.platform(chirp_end_time);
                obj.Attacker.platform(chirp_end_time);
    
                %repeat the process to obtain the positions and velocities at
                %the start of the given chirp
                [victim_pos, victim_vel] = obj.Victim.platform(chirp_end_time);
                [tgt_pos,tgt_vel] = obj.SimulatedTarget.platform(chirp_end_time);
                [attacker_pos, attacker_vel] = obj.Attacker.platform(chirp_end_time);
            end
        end     
        
        function received_signal = generate_sensing_subsystem_received_signal(obj)
            %{
                Purpose: the following code is used to determine the signal
                that the sensing subsystem would receive from the victim.
                Outputs:
                    received_signal: the signal received from the victim
            %}
            received_signal = zeros(obj.Attacker.Subsystem_spectrum_sensing.spectogram_params.num_ADC_samples_per_spectogram,1);
            received_signal_assembled = 0;
        
            while ~received_signal_assembled
                %if all of the chirps for a frame have been sent, just send
                %zeros until the end of the current frame
                if obj.sensing_subsystem_support.chirp_count >= obj.Victim.NumChirps
                    
                    if obj.sensing_subsystem_support.received_samples_left ...
                            <= obj.sensing_subsystem_support.frame_samples_left   %done sending chirps, but not end of the frame
                        points_to_insert = obj.sensing_subsystem_support.received_samples_left;
                        
        
                        %the following code optimizes the while loop script so that
                        %no computations are made when the "victim" is idle and not
                        %transmitting. To run the full thing, just comment out the
                        %if statement and only leave what is below the if
                        %statement.
                        if obj.sensing_subsystem_support.received_samples_left ...
                                ~= obj.sensing_subsystem_support.max_received_sample_index
                            received_signal(obj.sensing_subsystem_support.received_sample_index:...
                                obj.sensing_subsystem_support.received_sample_index + points_to_insert - 1)...
                                = zeros(points_to_insert,1);
                            received_signal_assembled = 1;
                        else
                            obj.Attacker.Subsystem_spectrum_sensing.sampling_window_count =...
                                obj.Attacker.Subsystem_spectrum_sensing.sampling_window_count + 1;
                            obj.Attacker.Subsystem_spectrum_sensing.previous_spectogram_points = [];
                        end
                        
                        obj.sensing_subsystem_support.received_sample_index = 1;
                        obj.sensing_subsystem_support.received_samples_left = obj.sensing_subsystem_support.max_received_sample_index;
        
                        obj.sensing_subsystem_support.frame_samples_left = ...
                            obj.sensing_subsystem_support.frame_samples_left -...
                            points_to_insert;
                        
                        if obj.sensing_subsystem_support.frame_samples_left == 0
                            obj.sensing_subsystem_support.frame_samples_left = ...
                                obj.sensing_subsystem_support.max_frame_sample_index;
                            obj.sensing_subsystem_support.frame_count = ...
                                obj.sensing_subsystem_support.frame_count + 1;
                            obj.sensing_subsystem_support.chirp_count = 0;
                        end
        
                    else                                                %represents the end of the specified frame
                        points_to_insert = ...
                            obj.sensing_subsystem_support.frame_samples_left;
                        
                        received_signal(obj.sensing_subsystem_support.received_sample_index:...
                            obj.sensing_subsystem_support.received_sample_index + points_to_insert - 1) = zeros(points_to_insert,1);
                        
                        
                        %increment the frame counter and start sending chirps again
                        %as new frame has started
                        obj.sensing_subsystem_support.frame_samples_left = ...
                            obj.sensing_subsystem_support.max_frame_sample_index;
                        obj.sensing_subsystem_support.frame_count = ...
                            obj.sensing_subsystem_support.frame_count + 1;
        
                        obj.sensing_subsystem_support.chirp_count = 0;
        
                        obj.sensing_subsystem_support.received_sample_index = ...
                            obj.sensing_subsystem_support.received_sample_index + points_to_insert;
                        obj.sensing_subsystem_support.received_samples_left ...
                            = obj.sensing_subsystem_support.max_received_sample_index -...
                            obj.sensing_subsystem_support.received_sample_index + 1;
                    end
                    
                %the next two statements deal with assembling chirps to be sent out

                %if there are more received samples to be added than chirp
                %samples to be added
                elseif obj.sensing_subsystem_support.received_samples_left >= ...
                        obj.sensing_subsystem_support.chirp_samples_left
                    points_to_insert = obj.sensing_subsystem_support.chirp_samples_left;
                    received_signal(obj.sensing_subsystem_support.received_sample_index:...
                        obj.sensing_subsystem_support.received_sample_index + points_to_insert - 1) =...
                        obj.Victim.chirp(obj.sensing_subsystem_support.chirp_sample_index:...
                        obj.sensing_subsystem_support.chirp_sample_index + points_to_insert - 1);
                    
                    obj.sensing_subsystem_support.chirp_sample_index = 1;
                    obj.sensing_subsystem_support.chirp_samples_left = ...
                        obj.sensing_subsystem_support.max_chirp_sample_index;
                    obj.sensing_subsystem_support.chirp_count = ...
                        obj.sensing_subsystem_support.chirp_count + 1;
                
                    obj.sensing_subsystem_support.received_sample_index = ...
                        obj.sensing_subsystem_support.received_sample_index + points_to_insert;
                    obj.sensing_subsystem_support.received_samples_left = ...
                        obj.sensing_subsystem_support.max_received_sample_index - ...
                        obj.sensing_subsystem_support.received_sample_index + 1;
                    
                    obj.sensing_subsystem_support.frame_samples_left = ...
                        obj.sensing_subsystem_support.frame_samples_left - points_to_insert;
        
        
                    if obj.sensing_subsystem_support.received_samples_left <= 0
                        obj.sensing_subsystem_support.received_sample_index = 1;
                        obj.sensing_subsystem_support.received_samples_left = ...
                            obj.sensing_subsystem_support.max_received_sample_index;
                        received_signal_assembled = 1;
                    end
                    
                %otherwise, there must be more chirp samples to insert for
                %the current chirp than there are remaining available
                %received samples
                else
                    points_to_insert = obj.sensing_subsystem_support.received_samples_left;
                    received_signal(obj.sensing_subsystem_support.received_sample_index:...
                        obj.sensing_subsystem_support.received_sample_index + points_to_insert - 1) =...
                        obj.Victim.chirp(obj.sensing_subsystem_support.chirp_sample_index:...
                        obj.sensing_subsystem_support.chirp_sample_index + points_to_insert - 1);
                
                    obj.sensing_subsystem_support.chirp_sample_index = ...
                        obj.sensing_subsystem_support.chirp_sample_index + points_to_insert;
                    obj.sensing_subsystem_support.chirp_samples_left = ...
                        obj.sensing_subsystem_support.max_chirp_sample_index - ...
                        obj.sensing_subsystem_support.chirp_sample_index + 1;
                
                    obj.sensing_subsystem_support.received_sample_index = 1;
                    obj.sensing_subsystem_support.received_samples_left = ...
                        obj.sensing_subsystem_support.max_received_sample_index;
                    received_signal_assembled = 1;
        
                    obj.sensing_subsystem_support.frame_samples_left =...
                        obj.sensing_subsystem_support.frame_samples_left - points_to_insert;
                end
            end
        end

    end

end