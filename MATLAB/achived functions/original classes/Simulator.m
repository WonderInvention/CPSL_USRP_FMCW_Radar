%% ARCHIVED, DO NOT USE EXCEPT FOR REFERENCE

classdef Simulator < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties (Access = public)
        Attacker
        Defender
        SimulatedTarget

        %variables used in computing the FMCW simulation - MOVE TO PRIVATE
        channel_target          %phased.FreeSpace object for a target
        channel_attacker        %phased.FreeSpace object for an attacker
    end

    properties (Access = private)

        %physical constants
        c

        
    end

    methods (Access = public)
        function obj = Simulator()
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            obj.Attacker = Radar();
            obj.Defender = Radar();
            obj.c = physconst('LightSpeed');
        end

%% [1] Functions to configure the simulation

        function load_default_defender_params(obj)
            %setup the defender's chirp parameters
            obj.Defender.StartFrequency_GHz         = 77.0;
            obj.Defender.FrequencySlope_MHz_us      = 65.0;
            obj.Defender.IdleTime_us                = 7.0;
            obj.Defender.TxStartTime_us             = 0;
            obj.Defender.ADC_ValidStartTime_us      = 6.3;
            obj.Defender.ADC_Samples                = 224;
            obj.Defender.ADC_SampleRate_MSps        = 4.1830;
            obj.Defender.RampEndTime_us             = 60.85;

            %setup the defender's frame parameters
            obj.Defender.NumChirps                  = 16;
            obj.Defender.FramePeriodicity_ms        = 33.33;
            
            %define plot color default values
            obj.Defender.plotResolution_us = 0.01;
            obj.Defender.tx_period_plot_color = 'blue';
            obj.Defender.tx_sampling_period_plot_color = 'cyan';
            obj.Defender.radar_name = 'Defender';

            %set timing offset to zero as this is the defender
            obj.Defender.timing_offset_us = 0;

            %compute all remaining "calculated" values
            obj.Defender.compute_calculated_vals();

            %set defender to "enabled"
            obj.Defender.enabled = true;
        end

        function load_gesture_defender_params(obj)
            %setup the defender's chirp parameters
            obj.Defender.StartFrequency_GHz         = 77.0;
            obj.Defender.FrequencySlope_MHz_us      = 26.67;
            obj.Defender.IdleTime_us                = 10;
            obj.Defender.TxStartTime_us             = 0;
            obj.Defender.ADC_ValidStartTime_us      = 6.3;
            obj.Defender.ADC_Samples                = 512;
            obj.Defender.ADC_SampleRate_MSps        = 3.718;
            obj.Defender.RampEndTime_us             = 150;
            
            %setup the defender's frame parameters
            obj.Defender.NumChirps                  = 255;
            obj.Defender.FramePeriodicity_ms        = 40.96;
            
            %define plot color default values
            obj.Defender.plotResolution_us = 0.01;
            obj.Defender.tx_period_plot_color = 'blue';
            obj.Defender.tx_sampling_period_plot_color = 'cyan';
            obj.Defender.radar_name = 'Defender';

            %set timing offset to zero as this is the defender
            obj.Defender.timing_offset_us = 0;

            %compute all remaining "calculated" values
            obj.Defender.compute_calculated_vals();

            %set defender to "enabled"
            obj.Defender.enabled = true;
        end

        function load_default_attacker_params(obj)
            %setup the defender's chirp parameters
            obj.Attacker.StartFrequency_GHz         = 77.0;
            obj.Attacker.FrequencySlope_MHz_us      = 65.0;
            obj.Attacker.IdleTime_us                = 7.0;
            obj.Attacker.TxStartTime_us             = 0;
            obj.Attacker.ADC_ValidStartTime_us      = 6.3;
            obj.Attacker.ADC_Samples                = 256;
            obj.Attacker.ADC_SampleRate_MSps        = 13.0445;
            obj.Attacker.RampEndTime_us             = 26.925;
            
            %setup the defender's frame parameters
            obj.Attacker.NumChirps                  = 32;
            obj.Attacker.FramePeriodicity_ms        = 33.4965;
            
            %define plot color default values
            obj.Attacker.plotResolution_us = 0.01;
            obj.Attacker.tx_period_plot_color = 'red';
            obj.Attacker.tx_sampling_period_plot_color = 'magenta';
            obj.Attacker.radar_name = 'Attacker';

            %define set the default offset to be 0us
            obj.Attacker.timing_offset_us = 0;

            %compute all remaining "calculated" values
            obj.Attacker.compute_calculated_vals();

            %set attacker to "enabled"
            obj.Attacker.enabled = true;
        end

        function load_attacker_scenario_2(obj)

            %setup the defender's chirp parameters
            obj.Attacker.StartFrequency_GHz         = 77.0;
            obj.Attacker.FrequencySlope_MHz_us      = 20.0;
            obj.Attacker.IdleTime_us                = 7.0;
            obj.Attacker.TxStartTime_us             = 0;
            obj.Attacker.ADC_ValidStartTime_us      = 6.3;
            obj.Attacker.ADC_Samples                = 256;
            obj.Attacker.ADC_SampleRate_MSps        = 13.0445;
            obj.Attacker.RampEndTime_us             = 29.625;
            
            %setup the defender's frame parameters
            obj.Attacker.NumChirps                  = 32;
            obj.Attacker.FramePeriodicity_ms        = 33.4965;
            
            %define plot color default values
            obj.Attacker.plotResolution_us = 0.01;
            obj.Attacker.tx_period_plot_color = 'red';
            obj.Attacker.tx_sampling_period_plot_color = 'magenta';
            obj.Attacker.radar_name = 'Attacker';

            %define set the default offset to be 0us
            obj.Attacker.timing_offset_us = 0;

            %compute all remaining "calculated" values
            obj.Attacker.compute_calculated_vals();

            %set attacker to "enabled"
            obj.Attacker.enabled = true;
        end

        function load_target_default(obj)
            %{
                Purpose: load a default target and set it to the "enabled"
                state
            %}

            %configure the simulated target
            distance_m = 5;
            speed_meters_per_s = 7;
            rcs_sq_meters = db2pow(min(10*log10(distance_m)+5,20));
            operating_frequency_Hz = obj.Defender.StartFrequency_GHz * 1e9;
            enabled = true;

            obj.SimulatedTarget = Target(distance_m,speed_meters_per_s,rcs_sq_meters,operating_frequency_Hz,enabled);
        end

        function load_default_attacker_and_defender_position_and_velocity(obj)
            %{
                Purpose: configures a default scenario for the attacker and
                defender positions and velocities
            %}
            obj.Defender.position_m = [0;0;0];
            obj.Defender.velocity_m_per_s = [0;0;0];
            obj.Defender.platform = phased.Platform( ...
                'InitialPosition',obj.Defender.position_m, ...
                'Velocity',obj.Defender.velocity_m_per_s);

            obj.Attacker.position_m = [7;0;0];
            obj.Attacker.velocity_m_per_s = [0;0;0];
            obj.Attacker.platform = phased.Platform( ...
                'InitialPosition',obj.Attacker.position_m, ...
                'Velocity',obj.Attacker.velocity_m_per_s);
        end
        
        function configure_FMCW_parameters(obj)
            %{
                Purpose: initializes all of the needed waveform parameters
                    to simulate the actual waveforms for the attacker and
                    defender
                TODO: code is currently only configured to plot the defender's
                    range fft's and whatnot. The attacker plots aren't generated
                    at this time, but would be helpful to have in the future
            %}

            obj.Defender.downsample_factor = ceil(2 * 4e9 / (obj.Defender.ADC_SampleRate_MSps * 1e6));

            obj.Defender.FMCW_sampling_rate_Hz = obj.Defender.ADC_SampleRate_MSps * 1e6 * obj.Defender.downsample_factor;
            obj.Attacker.FMCW_sampling_rate_Hz = obj.Defender.FMCW_sampling_rate_Hz;

            %The FMCWWaveform function only starts from 0Hz, if the start
            %frequency isn't 77.0 GHz, we need to determine the additional
            %BW and sweep time required. Then, when we generate the
            %waveform, we will cutoff the extra signal that we added.

                %for defender
                if obj.Defender.StartFrequency_GHz ~= 77.0
                    obj.Defender.additional_BW_MHz = (obj.Defender.StartFrequency_GHz - 77.0) * 1e3;
                    obj.Defender.additional_sweep_time_us = obj.Defender.additional_BW_MHz / obj.Defender.FrequencySlope_MHz_us;
                else
                    obj.Defender.additional_BW_MHz = 0;
                    obj.Defender.additional_sweep_time_us = 0;
                end

                %for attacker
                if obj.Attacker.StartFrequency_GHz ~= 77.0
                    obj.Attacker.additional_BW_MHz = (obj.Attacker.StartFrequency_GHz - 77.0) * 1e3;
                    obj.Attacker.additional_sweep_time_us = obj.Attacker.additional_BW_MHz / obj.Attacker.FrequencySlope_MHz_us;
                else
                    obj.Attacker.additional_BW_MHz = 0;
                    obj.Attacker.additional_sweep_time_us = 0;
                end

            %set sweep time
            obj.Defender.sweep_time = ceil(obj.Defender.FMCW_sampling_rate_Hz * (obj.Defender.RampEndTime_us + obj.Defender.additional_sweep_time_us) * 1e-6) / obj.Defender.FMCW_sampling_rate_Hz;
            obj.Attacker.sweep_time = ceil(obj.Attacker.FMCW_sampling_rate_Hz * (obj.Attacker.RampEndTime_us + obj.Attacker.additional_sweep_time_us) * 1e-6) / obj.Attacker.FMCW_sampling_rate_Hz;

            %configure waveforms
            obj.Defender.waveform = phased.FMCWWaveform( ...
                'SampleRate', obj.Defender.FMCW_sampling_rate_Hz, ...
                'SweepTime', obj.Defender.sweep_time,...
                'SweepBandwidth', (obj.Defender.Chirp_Tx_Bandwidth_MHz + obj.Defender.additional_BW_MHz) * 1e6, ...
                'SweepDirection', 'Up',...
                'SweepInterval', 'Positive',...
                'OutputFormat','Sweeps',...
                'NumSweeps',1);

            obj.Attacker.waveform = phased.FMCWWaveform( ...
                'SampleRate', obj.Attacker.FMCW_sampling_rate_Hz, ...
                'SweepTime', obj.Attacker.sweep_time,...
                'SweepBandwidth', (obj.Attacker.Chirp_Tx_Bandwidth_MHz + obj.Attacker.additional_BW_MHz) * 1e6, ...
                'SweepDirection', 'Up',...
                'SweepInterval', 'Positive',...
                'OutputFormat','Sweeps',...
                'NumSweeps',1);
            
            obj.channel_target = phased.FreeSpace( ...
                "PropagationSpeed",obj.c, ...
                "OperatingFrequency", obj.Defender.StartFrequency_GHz * 1e9, ...
                "SampleRate", obj.Defender.FMCW_sampling_rate_Hz, ...
                "TwoWayPropagation", true);

            obj.channel_attacker = phased.FreeSpace( ...
                "PropagationSpeed",obj.c, ...
                "OperatingFrequency", obj.Defender.StartFrequency_GHz * 1e9, ...
                "SampleRate", obj.Defender.FMCW_sampling_rate_Hz, ...
                "TwoWayPropagation", false);

            obj.Attacker.configure_transmitter_and_receiver();
            obj.Attacker.configure_lowpass_filter();

            obj.Defender.configure_transmitter_and_receiver();
            obj.Defender.configure_lowpass_filter();
        end

        %realistic attacker and victim scenario
        function load_realistic_defender_params(obj)
            %setup the defender's chirp parameters
            obj.Defender.StartFrequency_GHz         = 77.0;
            obj.Defender.FrequencySlope_MHz_us      = 10.76;
            obj.Defender.IdleTime_us                = 7.0;
            obj.Defender.TxStartTime_us             = 0;
            obj.Defender.ADC_ValidStartTime_us      = 6.3;
            obj.Defender.ADC_Samples                = 256;
            obj.Defender.ADC_SampleRate_MSps        = 7.17;
            obj.Defender.RampEndTime_us             = 43;

            %setup the defender's frame parameters
            obj.Defender.NumChirps                  = 32;
            obj.Defender.FramePeriodicity_ms        = 33.33;
            
            %define plot color default values
            obj.Defender.plotResolution_us = 0.01;
            obj.Defender.tx_period_plot_color = 'blue';
            obj.Defender.tx_sampling_period_plot_color = 'cyan';
            obj.Defender.radar_name = 'Defender';

            %set timing offset to zero as this is the defender
            obj.Defender.timing_offset_us = 0;

            %compute all remaining "calculated" values
            obj.Defender.compute_calculated_vals();

            %set defender to "enabled"
            obj.Defender.enabled = true;
        end

        function load_realisitc_attacker_params(obj)
            %setup the defender's chirp parameters
            obj.Attacker.StartFrequency_GHz         = 77.0;
            obj.Attacker.FrequencySlope_MHz_us      = 10.76;
            obj.Attacker.IdleTime_us                = 7.0;
            obj.Attacker.TxStartTime_us             = 0;
            obj.Attacker.ADC_ValidStartTime_us      = 6.3;
            obj.Attacker.ADC_Samples                = 256;
            obj.Attacker.ADC_SampleRate_MSps        = 7.17;
            obj.Attacker.RampEndTime_us             = 43;
            
            %setup the defender's frame parameters
            obj.Attacker.NumChirps                  = 32;
            obj.Attacker.FramePeriodicity_ms        = 33.4965;
            
            %define plot color default values
            obj.Attacker.plotResolution_us = 0.01;
            obj.Attacker.tx_period_plot_color = 'red';
            obj.Attacker.tx_sampling_period_plot_color = 'magenta';
            obj.Attacker.radar_name = 'Attacker';

            %define set the default offset to be 0us
            obj.Attacker.timing_offset_us = 0;

            %compute all remaining "calculated" values
            obj.Attacker.compute_calculated_vals();

            %set attacker to "enabled"
            obj.Attacker.enabled = true;
        end

        function load_target_realistic(obj)
            %{
                Purpose: load a default target and set it to the "enabled"
                state
            %}

            %configure the simulated target
            distance_m = 55;
            speed_meters_per_s = 7;
            rcs_sq_meters = db2pow(min(10*log10(distance_m)+5,20));
            operating_frequency_Hz = obj.Defender.StartFrequency_GHz * 1e9;
            enabled = true;

            obj.SimulatedTarget = Target(distance_m,speed_meters_per_s,rcs_sq_meters,operating_frequency_Hz,enabled);
        end
        
        function load_realistic_attacker_and_defender_position_and_velocity(obj)
            %{
                Purpose: configures a default scenario for the attacker and
                defender positions and velocities
            %}
            obj.Defender.position_m = [0;0;0];
            obj.Defender.velocity_m_per_s = [0;0;0];
            obj.Defender.platform = phased.Platform( ...
                'InitialPosition',obj.Defender.position_m, ...
                'Velocity',obj.Defender.velocity_m_per_s);

            obj.Attacker.position_m = [75;0;0];
            obj.Attacker.velocity_m_per_s = [0;0;0];
            obj.Attacker.platform = phased.Platform( ...
                'InitialPosition',obj.Attacker.position_m, ...
                'Velocity',obj.Attacker.velocity_m_per_s);
        end

        %low bandwidth attacks
            %for realistic defender
                function load_realisitc_attacker_params_40MHz(obj)
                    %setup the defender's chirp parameters
                    obj.Attacker.StartFrequency_GHz         = 77.1076;
                    obj.Attacker.FrequencySlope_MHz_us      = 10.76;
                    obj.Attacker.IdleTime_us                = 46.28;
                    obj.Attacker.TxStartTime_us             = 0;
                    obj.Attacker.ADC_ValidStartTime_us      = 0;
                    obj.Attacker.ADC_Samples                = 26;
                    obj.Attacker.ADC_SampleRate_MSps        = 7.17;
                    obj.Attacker.RampEndTime_us             = 3.72;
                    
                    %setup the defender's frame parameters
                    obj.Attacker.NumChirps                  = 32;
                    obj.Attacker.FramePeriodicity_ms        = 33.33;
                    
                    %define plot color default values
                    obj.Attacker.plotResolution_us = 0.01;
                    obj.Attacker.tx_period_plot_color = 'red';
                    obj.Attacker.tx_sampling_period_plot_color = 'magenta';
                    obj.Attacker.radar_name = 'Attacker';
        
                    %define set the default offset to be 0us
                    obj.Attacker.timing_offset_us = -29.28;
        
                    %compute all remaining "calculated" values
                    obj.Attacker.compute_calculated_vals();
        
                    %set attacker to "enabled"
                    obj.Attacker.enabled = true;
                end

                function load_realisitc_attacker_params_80MHz(obj)
                    %setup the defender's chirp parameters
                    obj.Attacker.StartFrequency_GHz         = 77.1076;
                    obj.Attacker.FrequencySlope_MHz_us      = 10.76;
                    obj.Attacker.IdleTime_us                = 42.57;
                    obj.Attacker.TxStartTime_us             = 0;
                    obj.Attacker.ADC_ValidStartTime_us      = 0;
                    obj.Attacker.ADC_Samples                = 53;
                    obj.Attacker.ADC_SampleRate_MSps        = 7.17;
                    obj.Attacker.RampEndTime_us             = 7.43;
                    
                    %setup the defender's frame parameters
                    obj.Attacker.NumChirps                  = 32;
                    obj.Attacker.FramePeriodicity_ms        = 33.33;
                    
                    %define plot color default values
                    obj.Attacker.plotResolution_us = 0.01;
                    obj.Attacker.tx_period_plot_color = 'red';
                    obj.Attacker.tx_sampling_period_plot_color = 'magenta';
                    obj.Attacker.radar_name = 'Attacker';
        
                    %define set the default offset to be 0us
                    obj.Attacker.timing_offset_us = -25.57;
        
                    %compute all remaining "calculated" values
                    obj.Attacker.compute_calculated_vals();
        
                    %set attacker to "enabled"
                    obj.Attacker.enabled = true;
                end

                function load_realisitc_attacker_params_150MHz(obj)
                    %setup the defender's chirp parameters
                    obj.Attacker.StartFrequency_GHz         = 77.1076;
                    obj.Attacker.FrequencySlope_MHz_us      = 10.76;
                    obj.Attacker.IdleTime_us                = 36.06;
                    obj.Attacker.TxStartTime_us             = 0;
                    obj.Attacker.ADC_ValidStartTime_us      = 0;
                    obj.Attacker.ADC_Samples                = 99;
                    obj.Attacker.ADC_SampleRate_MSps        = 7.17;
                    obj.Attacker.RampEndTime_us             = 13.94;
                    
                    %setup the defender's frame parameters
                    obj.Attacker.NumChirps                  = 32;
                    obj.Attacker.FramePeriodicity_ms        = 33.33;
                    
                    %define plot color default values
                    obj.Attacker.plotResolution_us = 0.01;
                    obj.Attacker.tx_period_plot_color = 'red';
                    obj.Attacker.tx_sampling_period_plot_color = 'magenta';
                    obj.Attacker.radar_name = 'Attacker';
        
                    %define set the default offset to be 0us
                    obj.Attacker.timing_offset_us = -19.06;
        
                    %compute all remaining "calculated" values
                    obj.Attacker.compute_calculated_vals();
        
                    %set attacker to "enabled"
                    obj.Attacker.enabled = true;
                end
            
                %for default defender (~4 GHz BW)
                function load_default_attacker_params_80MHz(obj)
                    %setup the defender's chirp parameters
                    obj.Attacker.StartFrequency_GHz         = 77.65;
                    obj.Attacker.FrequencySlope_MHz_us      = 65.0;
                    obj.Attacker.IdleTime_us                = 66.6;
                    obj.Attacker.TxStartTime_us             = 0;
                    obj.Attacker.ADC_ValidStartTime_us      = 0;
                    obj.Attacker.ADC_Samples                = 5;
                    obj.Attacker.ADC_SampleRate_MSps        = 4.183;
                    obj.Attacker.RampEndTime_us             = 1.25;
                    
                    %setup the defender's frame parameters
                    obj.Attacker.NumChirps                  = 32;
                    obj.Attacker.FramePeriodicity_ms        = 33.33;
                    
                    %define plot color default values
                    obj.Attacker.plotResolution_us = 0.01;
                    obj.Attacker.tx_period_plot_color = 'red';
                    obj.Attacker.tx_sampling_period_plot_color = 'magenta';
                    obj.Attacker.radar_name = 'Attacker';
        
                    %define set the default offset to be 0us
                    obj.Attacker.timing_offset_us = -49.6;
        
                    %compute all remaining "calculated" values
                    obj.Attacker.compute_calculated_vals();
        
                    %set attacker to "enabled"
                    obj.Attacker.enabled = true;
                end

                function load_default_attacker_params_150MHz(obj)
                    %setup the defender's chirp parameters
                    obj.Attacker.StartFrequency_GHz         = 77.65;
                    obj.Attacker.FrequencySlope_MHz_us      = 65.0;
                    obj.Attacker.IdleTime_us                = 65.54;
                    obj.Attacker.TxStartTime_us             = 0;
                    obj.Attacker.ADC_ValidStartTime_us      = 0;
                    obj.Attacker.ADC_Samples                = 9;
                    obj.Attacker.ADC_SampleRate_MSps        = 4.183;
                    obj.Attacker.RampEndTime_us             = 2.31;
                    
                    %setup the defender's frame parameters
                    obj.Attacker.NumChirps                  = 32;
                    obj.Attacker.FramePeriodicity_ms        = 33.33;
                    
                    %define plot color default values
                    obj.Attacker.plotResolution_us = 0.01;
                    obj.Attacker.tx_period_plot_color = 'red';
                    obj.Attacker.tx_sampling_period_plot_color = 'magenta';
                    obj.Attacker.radar_name = 'Attacker';
        
                    %define set the default offset to be 0us
                    obj.Attacker.timing_offset_us = -48.54;
        
                    %compute all remaining "calculated" values
                    obj.Attacker.compute_calculated_vals();
        
                    %set attacker to "enabled"
                    obj.Attacker.enabled = true;
                end

                function load_default_attacker_params_400MHz(obj)
                    %setup the defender's chirp parameters
                    obj.Attacker.StartFrequency_GHz         = 77.65;
                    obj.Attacker.FrequencySlope_MHz_us      = 65.0;
                    obj.Attacker.IdleTime_us                = 61.7;
                    obj.Attacker.TxStartTime_us             = 0;
                    obj.Attacker.ADC_ValidStartTime_us      = 0;
                    obj.Attacker.ADC_Samples                = 25;
                    obj.Attacker.ADC_SampleRate_MSps        = 4.183;
                    obj.Attacker.RampEndTime_us             = 6.15;
                    
                    %setup the defender's frame parameters
                    obj.Attacker.NumChirps                  = 32;
                    obj.Attacker.FramePeriodicity_ms        = 33.33;
                    
                    %define plot color default values
                    obj.Attacker.plotResolution_us = 0.01;
                    obj.Attacker.tx_period_plot_color = 'red';
                    obj.Attacker.tx_sampling_period_plot_color = 'magenta';
                    obj.Attacker.radar_name = 'Attacker';
        
                    %define set the default offset to be 0us
                    obj.Attacker.timing_offset_us = -44.7;
        
                    %compute all remaining "calculated" values
                    obj.Attacker.compute_calculated_vals();
        
                    %set attacker to "enabled"
                    obj.Attacker.enabled = true;
                end

%% [2] Functions to generate frequency over time for a given chirp or frame

        function [radar_tx,radar_f_tx, radar_t_sampling, radar_f_sampling] = plot_specified_radar_chirp(obj, axis, radar, radar_frame, radar_chirp)
            %{
                Purpose: plot a desired chirp in a given frame for a
                    specified radar object
                Inputs:
                    axis: the axis object corresponding to the figure where
                    the plot will go
                    radar: the radar object whose chirp is to be plotted
                    radar_frame: the desired frame to plot
                    radar_chirp: the desired chirp to plot
                Outputs:
                    [radar_tx,radar_f_tx]: the [t,f] values corresponding
                        to the desired chirp tx period
                    [radar_t_sampling, radar_f_sampling]: the [t,f] values
                        corresponding to the the desired chirp's sampling
                        period
            %}
            frame_start_us = radar.FramePeriodicity_ms * 1000 * (radar_frame - 1) + radar.timing_offset_us;
            chirp_start_us = frame_start_us + radar.ChirpCycleTime_us * (radar_chirp - 1);
            chirp_end_us = frame_start_us + radar.ChirpCycleTime_us * (radar_chirp);
            
            [radar_tx,radar_f_tx] = radar.generate_chirp_f_over_t_vals(chirp_start_us);
            [radar_t_sampling,radar_f_sampling] = radar.generate_chirp_f_over_Tsampling_vals(chirp_start_us);

            %plot the desired chirp
            obj.plot_radar_chirp(axis,radar,radar_tx,radar_f_tx,radar_t_sampling,radar_f_sampling);
        end

        function plot_attacker_and_defender_chirps(obj,axis,defender_frame,defender_chirp)
            %{
                Purpose: generate a plot for a specific defender chirp in a
                    given frame and overlay any attacker chirps that would
                    overlap with that defender chirp
                Inputs:
                    axis: the axis object corresponding to the figure to
                        plot in
                    defender_frame: the desired frame number of the
                        defender
                    defender_chirp: the desired chirp number
            %}

            %calculate the [t,f] values for the defender chirp

            defender_frame_start_us = obj.Defender.FramePeriodicity_ms * 1000 * (defender_frame - 1) + obj.Defender.timing_offset_us;
            defender_chirp_start_us = defender_frame_start_us + obj.Defender.ChirpCycleTime_us * (defender_chirp - 1);
            defender_chirp_end_us = defender_frame_start_us + obj.Defender.ChirpCycleTime_us * (defender_chirp);
            
            [defender_tx,defender_f_tx] = obj.Defender.generate_chirp_f_over_t_vals(defender_chirp_start_us);
            [defender_t_sampling,defender_f_sampling] = obj.Defender.generate_chirp_f_over_Tsampling_vals(defender_chirp_start_us);

            %now, calculate the [t,f] values for the attacker
            
            %determine the frames that will need to be plotted
            attacker_frame = floor((defender_chirp_start_us - obj.Attacker.timing_offset_us)/(obj.Attacker.FramePeriodicity_ms*1000))+1;
            attacker_frame_start_us = obj.Attacker.FramePeriodicity_ms * 1000 * (attacker_frame - 1) + obj.Attacker.timing_offset_us;


            %initialize the attacker chirp to start with
            attacker_chirp = floor((defender_chirp_start_us - attacker_frame_start_us)/ obj.Attacker.ChirpCycleTime_us)+1;
            if attacker_chirp > obj.Attacker.NumChirps
                attacker_chirp = obj.Attacker.NumChirps;
            end

            %calculate the [t,f] values for the first chirp
            attacker_chirp_start_us = attacker_frame_start_us + obj.Attacker.ChirpCycleTime_us * (attacker_chirp - 1);
            attacker_chirp_end_us = attacker_frame_start_us + obj.Attacker.ChirpCycleTime_us * (attacker_chirp);
            [attacker_tx,attacker_f_tx] = obj.Attacker.generate_chirp_f_over_t_vals(attacker_chirp_start_us);
            [attacker_t_sampling,attacker_f_sampling] = obj.Attacker.generate_chirp_f_over_Tsampling_vals(attacker_chirp_start_us);

            while attacker_chirp_end_us < defender_chirp_end_us
                %move to the next chirp, but if at the last chirp in the frame, move to the first chirp
                %of the next frame. 
                if attacker_chirp < obj.Attacker.NumChirps
                    attacker_chirp = attacker_chirp + 1;
                else
                    attacker_frame = attacker_frame + 1;
                    attacker_frame_start_us = obj.Attacker.FramePeriodicity_ms * 1000 * (attacker_frame - 1) + obj.Attacker.timing_offset_us;
                    attacker_chirp = 1;
                end
                attacker_chirp_start_us = attacker_frame_start_us + obj.Attacker.ChirpCycleTime_us * (attacker_chirp - 1);
                attacker_chirp_end_us = attacker_frame_start_us + obj.Attacker.ChirpCycleTime_us * (attacker_chirp);

                %calculate the [t,f] values for this additional chirp
                [tx,f_tx] = obj.Attacker.generate_chirp_f_over_t_vals(attacker_chirp_start_us);
                [t_sampling,f_sampling] = obj.Attacker.generate_chirp_f_over_Tsampling_vals(attacker_chirp_start_us);

                %add to the existing [t,f] values
                attacker_tx = [attacker_tx, tx];
                attacker_f_tx = [attacker_f_tx,f_tx];
                attacker_t_sampling = [attacker_t_sampling,t_sampling];
                attacker_f_sampling = [attacker_f_sampling,f_sampling];
            end

            %identify the points where the defender would actually observe the
            
           [attack_region_t, attack_region_f] = obj.identify_valid_chirp_interference(defender_t_sampling, defender_f_sampling, attacker_tx, attacker_f_tx);

            %plot everything
            obj.plot_radar_chirp(axis,obj.Defender,defender_tx,defender_f_tx,defender_t_sampling,defender_f_sampling);
            hold on
            obj.plot_radar_chirp(axis,obj.Attacker,attacker_tx,attacker_f_tx,attacker_t_sampling,attacker_f_sampling);
            
            hold on
            p = scatter(axis,attack_region_t,attack_region_f,'.');
            p.DisplayName = sprintf('Valid attack region ');
            p.MarkerEdgeColor = 'green';
            hold off
            
            % by default, zoom the axis to only focus on the defender chirp
            set(axis,'XLim',[defender_chirp_start_us - 10,defender_chirp_end_us + 10]);
        end
        
        function [radar_tx,radar_f_tx, radar_t_sampling, radar_f_sampling] = plot_specified_radar_frame(obj, axis, radar, radar_frame)
            %{
                Purpose: plot a desired framr for a
                    specified radar object
                Inputs:
                    axis: the axis object corresponding to the figure where
                    the plot will go
                    radar: the radar object whose chirp is to be plotted
                    radar_frame: the desired frame to plot
                Outputs:
                    [radar_tx,radar_f_tx]: the [t,f] values corresponding
                        to the desired frame tx period
                    [radar_t_sampling, radar_f_sampling]: the [t,f] values
                        corresponding to the the desired frame's sampling
                        period
            %}
            radar.plotResolution_us = 10 * radar.plotResolution_us;

            frame_start_ms = radar.FramePeriodicity_ms * (radar_frame - 1) + radar.timing_offset_us / 1000;
            
            [radar_tx,radar_f_tx] = radar.generate_frame_f_over_T_vals(frame_start_ms);

            [radar_t_sampling,radar_f_sampling] = radar.generate_frame_f_over_Tsampling_vals(frame_start_ms);

            %plot the desired chirp
            obj.plot_radar_frame(axis,radar,radar_tx,radar_f_tx,radar_t_sampling,radar_f_sampling);

            radar.plotResolution_us = radar.plotResolution_us / 10;
        end

        function plot_attacker_and_defender_frames(obj,axis,defender_frame)
                %{
                Purpose: generate a plot for a specific defender frame and 
                    overlay any attacker frames that would overlap with 
                    that defender chirp
                Inputs:
                    axis: the axis object corresponding to the figure to
                        plot in
                    defender_frame: the desired frame number of the
                        defender
                %}
             %calculate the [t,f] values for the defender chirp

            defender_frame_start_ms = obj.Defender.FramePeriodicity_ms * (defender_frame - 1) + obj.Defender.timing_offset_us / 1000;
            defender_frame_end_ms = defender_frame_start_ms + obj.Defender.FramePeriodicity_ms * (defender_frame);
            
            [defender_tx,defender_f_tx] = obj.Defender.generate_frame_f_over_T_vals(defender_frame_start_ms);
            [defender_t_sampling,defender_f_sampling] = obj.Defender.generate_frame_f_over_Tsampling_vals(defender_frame_start_ms);

            %now, calculate the [t,f] values for the attacker
            
            %determine the frames that will need to be plotted
            attacker_frame = floor((defender_frame_start_ms - (obj.Attacker.timing_offset_us/1000))/obj.Attacker.FramePeriodicity_ms) + 1;
            attacker_frame_start_ms = obj.Attacker.FramePeriodicity_ms * (attacker_frame - 1) + obj.Attacker.timing_offset_us / 1000;

            %calculate the [t,f] values for the first attacker frame
            attacker_frame_end_ms = attacker_frame_start_ms + obj.Attacker.FramePeriodicity_ms * (attacker_frame);
            [attacker_tx,attacker_f_tx] = obj.Attacker.generate_frame_f_over_T_vals(attacker_frame_start_ms);
            [attacker_t_sampling,attacker_f_sampling] = obj.Attacker.generate_frame_f_over_Tsampling_vals(attacker_frame_start_ms);

            while attacker_frame_end_ms < defender_frame_end_ms
                % if the attacker frame ends before the defender frame, 
                % move to the next attacker frame
                attacker_frame = attacker_frame + 1;
                attacker_frame_start_ms = obj.Attacker.FramePeriodicity_ms * (attacker_frame - 1) + obj.Attacker.timing_offset_us / 1000;
                attacker_frame_end_ms =  attacker_frame_start_ms + obj.Attacker.FramePeriodicity_ms * (attacker_frame);

                %calculate the [t,f] values for this additional frame
                [tx,f_tx] = obj.Attacker.generate_frame_f_over_T_vals(attacker_frame_start_ms);
                [t_sampling,f_sampling] = obj.Attacker.generate_frame_f_over_Tsampling_vals(attacker_frame_start_ms);

                %add to the existing [t,f] values
                attacker_tx = [attacker_tx, tx];
                attacker_f_tx = [attacker_f_tx,f_tx];
                attacker_t_sampling = [attacker_t_sampling,t_sampling];
                attacker_f_sampling = [attacker_f_sampling,f_sampling];
            end

            %identify the points where the defender would actually observe
            %the intersection of the attacker and defender frames
            
           [attack_region_t, attack_region_f] = obj.identify_valid_frame_interference(defender_t_sampling, defender_f_sampling, attacker_tx, attacker_f_tx);

            %plot everything
            obj.plot_radar_frame(axis,obj.Defender,defender_tx,defender_f_tx,defender_t_sampling,defender_f_sampling);
            hold on
            obj.plot_radar_frame(axis,obj.Attacker,attacker_tx,attacker_f_tx,attacker_t_sampling,attacker_f_sampling);
            
            hold on
            p = scatter(axis,attack_region_t,attack_region_f,'.');
            p.DisplayName = sprintf('Valid attack region ');
            p.MarkerEdgeColor = 'green';
            hold off

            % by default, zoom the axis to only focus on the defender frame
            set(axis,'XLim',[defender_frame_start_ms,defender_frame_start_ms + (obj.Defender.ChirpCycleTime_us * 1e-3 * obj.Defender.NumChirps)]);
        end

%% [3] FMCW simulations
%% [3.1] FMCW simulations - generate spectograms

        function FMCW_plot_radar_chirp_spectogram(obj, axis, radar, radar_frame, radar_chirp)
            
            %generate the signal
            frame_start_us = radar.FramePeriodicity_ms * 1000 * (radar_frame - 1) + radar.timing_offset_us;
            chirp_start_us = frame_start_us + radar.ChirpCycleTime_us * (radar_chirp - 1);
            
            [radar_t,radar_sig] = radar.FMCW_generate_chirp_sig_vals(chirp_start_us);

            %plot on the spectogram
            obj.plot_RF_spectogram(axis,radar_sig, radar.FMCW_sampling_rate_Hz, radar_t);

        end

        function plot_RF_spectogram(obj,axis,signal,sample_rate_Hz,t_chirp)
            %{
                Purpose: plot a spectogram plot on the given axis object
                Inputs:
                    axis: the axis to plot on
                    signal: the signal to plot for the spectogram
                    sample_rate_Hz: the sampling rate of the signal in Hz
                    t_chirp: the t-values corresponding to the chirp
            %}

            %define key spectogram parameters
            windowlength = 1028;
            noverlap = 512;
            nfft = 1028;

%             %plot the spectogram - previous method
%             axes(axis)
%             subplot(2,1,1);
%             spectrogram(signal,windowlength,noverlap,nfft,sample_rate_Hz, 'yaxis');
%             ylim([0 4.5])

            %plot the spectogram with correct time increments
            [s,f,t,ps] = spectrogram(signal,windowlength,noverlap,nfft,sample_rate_Hz);
            imagesc(axis,t_chirp,f * 1e-9 + 77.0, 10*log10(abs(ps)));
            title_str = sprintf('FMCW RF Waveform Spectogram');
            title(title_str);
            xlabel('Time (us)')
            ylabel('Frequency (GHz)')
            colorbar
            set(gca,'YDir','normal')
            set(gca,'YLim',[77,81.5])
        end

        function plot_IF_spectogram(obj,axis,signal,sample_rate_Hz)
            %{
                Purpose: plot a spectogram plot on the given axis object
                Inputs:
                    axis: the axis to plot on
                    signal: the signal to plot for the spectogram
                    sample_rate_Hz: the sampling rate of the signal in Hz
            %}

            %define key spectogram parameters
            windowlength = 64;
            noverlap = 32;
            nfft = 64;

            %plot the spectogram
            axes(axis)
            spectrogram(signal,windowlength,noverlap,nfft,sample_rate_Hz, 'yaxis');
            %ylim([0 20]);
        end
    
%% [3.2] FMCW simulations - support functions to simulate scene for a specific chirp

        function [defender_pos, defender_vel,attacker_pos, attacker_vel, tgt_pos,tgt_vel] = FMCW_determine_positions_and_velocities(obj,defender_frame,defender_chirp)
        %{
            Purpose: determine the position and velocity for the
                attacker,defender, and target at the start of a given chirp
                in a given frame
            Inputs:
                defender_chirp: the desired chirp
                defender_frame: the desired frame
            Outputs:
                [defender_pos, defender_vel]: the position and velocity of
                    the defender at the start of the desired chirp
                [attacker_pos, attacker_vel]: the position and velocity of
                    the attacker at the start of the desired chirp
                [target_pos, target_vel]: the position and velocity of
                    the attacker at the start of the desired chirp
        %}
            %calculate the time that the desired chirp will start at
            frame_start_time_s = (obj.Defender.FramePeriodicity_ms * 1e-3) * (defender_frame - 1);
            chirp_start_time_s = frame_start_time_s + (obj.Defender.ChirpCycleTime_us * 1e-6) * (defender_chirp -1);
            
            %reset the positions of each platform to be safe
            reset(obj.Defender.platform);
            reset(obj.Attacker.platform);
            reset(obj.SimulatedTarget.platform);

            if chirp_start_time_s <= 0
                time_increment = obj.Defender.ChirpCycleTime_us * 1e-6;
                [defender_pos, defender_vel] = obj.Defender.platform(time_increment);
                [tgt_pos,tgt_vel] = obj.SimulatedTarget.platform(time_increment);
                [attacker_pos, attacker_vel] = obj.Attacker.platform(time_increment);
            else
                %take a step for each platform so that they update to the
                %values at the start of the frame
                obj.Defender.platform(chirp_start_time_s);
                obj.SimulatedTarget.platform(chirp_start_time_s);
                obj.Attacker.platform(chirp_start_time_s);
    
                %repeat the process to obtain the positions and velocities at
                %the start of the given chirp
                [defender_pos, defender_vel] = obj.Defender.platform(chirp_start_time_s);
                [tgt_pos,tgt_vel] = obj.SimulatedTarget.platform(chirp_start_time_s);
                [attacker_pos, attacker_vel] = obj.Attacker.platform(chirp_start_time_s);
            end
        end
        
        function [defender_t,defender_sig, defender_chirp_start_us,defender_chirp_end_us] = FMCW_generate_defender_chirp_waveform(obj,defender_frame,defender_chirp)
            %{
                Purpose: simulate a specific defender chirp in a
                    given frame. Use the actual FMCW
                    waveform to simulate though
                Inputs:
                    defender_frame: the desired frame number of the
                        defender
                    defender_chirp: the desired chirp number
                Outputs:
                    [defender_t,defender_sig]: the [t,sig] values for the
                        defender's FMCW waveform
                    [defender_chirp_start_us,defender_chirp_end_us]: the
                        start and end time of the chirp in us
            %}
            %calculate the times and values for the defender chirp

            defender_frame_start_us = obj.Defender.FramePeriodicity_ms * 1000 * (defender_frame - 1) + obj.Defender.timing_offset_us;
            defender_chirp_start_us = defender_frame_start_us + obj.Defender.ChirpCycleTime_us * (defender_chirp - 1);
            defender_chirp_end_us = defender_frame_start_us + obj.Defender.ChirpCycleTime_us * (defender_chirp);
            
            [defender_t,defender_sig] = obj.Defender.FMCW_generate_chirp_sig_vals(defender_chirp_start_us);
        end

        function [attacker_t,attacker_sig] = FMCW_generate_attacker_chirp_waveforms(obj,defender_t,defender_sig, defender_chirp_start_us,defender_chirp_end_us)
            %{
                Purpose: Simulate an attacker chirps that would
                    overlap with that the given defender chirp. Use the 
                    actual FMCW waveform to simulate though.
                Inputs:
                    [defender_t,defender_sig]: the [t,sig] values for the
                        defender's FMCW waveform
                    [defender_chirp_start_us,defender_chirp_end_us]: the
                        start and end time of the chirp in us
                Outputs:
                    [attacker_t,attacker_sig]: the [t,sig] values for FMCW
                        waveforms corresponding from the attacker that
                        overlap with the given defender chirp
            %}

            %calculate the [t,sig] values for the attacker
             
            %determine the first frame that could potentially interfere
            attacker_frame = floor((defender_chirp_start_us - obj.Attacker.timing_offset_us)/(obj.Attacker.FramePeriodicity_ms*1000))+1;
            attacker_frame_start_us = obj.Attacker.FramePeriodicity_ms * 1000 * (attacker_frame - 1) + obj.Attacker.timing_offset_us;


            %initialize the attacker chirp to start with
            attacker_chirp = floor((defender_chirp_start_us - attacker_frame_start_us)/ obj.Attacker.ChirpCycleTime_us)+1;
            if attacker_chirp > obj.Attacker.NumChirps
                attacker_chirp = obj.Attacker.NumChirps;
            end

            %calculate the [t,sig] values for the first chirp
            attacker_chirp_start_us = attacker_frame_start_us + obj.Attacker.ChirpCycleTime_us * (attacker_chirp - 1);
            attacker_chirp_end_us = attacker_frame_start_us + obj.Attacker.ChirpCycleTime_us * (attacker_chirp);
            [attacker_t,attacker_sig] = obj.Attacker.FMCW_generate_chirp_sig_vals(attacker_chirp_start_us);

            [attacker_t,attacker_sig] = obj.FMCW_identify_valid_chirp_intersection(defender_t, defender_sig, attacker_t, attacker_sig);
            
            %determine how many zeros we need to fill in for the time
            %slots between the start of the victim chirp and the start of
            %the attacking chirp
            if size(attacker_t,1) >= 1 && attacker_t(1) - defender_t(1) > 0 %if there are points that overlap and the attacker starts after the defender
                increment = 1/obj.Defender.FMCW_sampling_rate_Hz * 1e6;
                t_additional = defender_t(1) : increment : attacker_t(1) - increment;
                t_additional = reshape(t_additional, size(t_additional,2),1);
                sig_additional = zeros(size(t_additional,1),1);

                %add to the existing [t,f] values
                attacker_t = [t_additional; attacker_t];
                attacker_sig = [sig_additional; attacker_sig];
            end

            %generate more attack chirps if necessary
            while attacker_chirp_end_us < defender_chirp_end_us
                %move to the next chirp, but if at the last chirp in the frame, move to the first chirp
                %of the next frame. 
                if attacker_chirp < obj.Attacker.NumChirps
                    attacker_chirp = attacker_chirp + 1;
                else
                    attacker_frame = attacker_frame + 1;
                    attacker_frame_start_us = obj.Attacker.FramePeriodicity_ms * 1000 * (attacker_frame - 1) + obj.Attacker.timing_offset_us;
                    attacker_chirp = 1;
                end
                attacker_chirp_start_us = attacker_frame_start_us + obj.Attacker.ChirpCycleTime_us * (attacker_chirp - 1);
                attacker_chirp_end_us = attacker_frame_start_us + obj.Attacker.ChirpCycleTime_us * (attacker_chirp);

                %calculate the [t,sig] values for this additional chirp
                [t,sig] = obj.Attacker.FMCW_generate_chirp_sig_vals(attacker_chirp_start_us);
                
                %determine the parts of the additional chirp that overlap
                %with the defender chirp
                [t,sig] = obj.FMCW_identify_valid_chirp_intersection(defender_t, defender_sig, t, sig);
                
                %determine how many zeros we need to fill in for the time
                %slots between chirps
                if size(t,1) >= 1 %if there are points to add to the attacker
                    increment = 1/obj.Defender.FMCW_sampling_rate_Hz * 1e6;
                    if size(attacker_t,1) >=1 %if there are already points for the attacker
                        t_additional = attacker_t(end) + increment : increment : t(1) - increment;
                        t_additional = reshape(t_additional, size(t_additional,2),1);
                        sig_additional = zeros(size(t_additional,1),1);
                    else % if there are no points for the attacker
                        t_additional = defender_t(1) : increment : t(1) - increment;
                        t_additional = reshape(t_additional, size(t_additional,2),1);
                        sig_additional = zeros(size(t_additional,1),1);
                    end

                    %add to the existing [t,f] values
                    attacker_t = [attacker_t;t_additional; t];
                    attacker_sig = [attacker_sig;sig_additional; sig];
                end

                
            end

            %fill in the rest of attacker instance with zeros so that its
            %the same length at the defender chirp
            increment = 1/obj.Defender.FMCW_sampling_rate_Hz * 1e6;
            if size(attacker_t,1) >= 1 %if the attacker had overlap with the defender
                t_additional = attacker_t(end) + increment : increment : defender_t(end);
                t_additional = reshape(t_additional, size(t_additional,2),1);
                sig_additional = zeros(size(t_additional,1),1);
            else % if the attacker did not have overlap with the defender
                t_additional = defender_t;
                sig_additional = zeros(size(defender_t,1),1);
            end

            %add to the existing [t,f] values
            attacker_t = [attacker_t;t_additional];
            attacker_sig = [attacker_sig;sig_additional];
        end
    
        function received_target_sig = FMCW_simulate_received_chirp_from_target (obj,defender_pos, defender_vel, tgt_pos,tgt_vel,defender_sig)
            %{
                Purpose: Simulate the signal that the defender would
                    receive from its transmit signal being reflected off 
                    of the target
                Inputs:
                    [defender_pos, defender_vel]: the position and velocity
                        vectors for the defender (generated from 
                        phased.platform)
                    [tgt_pos,tgt_vel]: the position and velocity vectors
                        for the target (generated from phased.platform)
                    defender_sig: the sig values of the
                        defender's FMCW chirp waveform
                Outputs:
                    target_sig: the FMCW waveform signal after it has been
                        received by the defender's receiver, but before
                        its been dechirped
            %}
            
            %transmit the defender chirp
            received_target_sig = obj.Defender.transmitter(defender_sig);

            %propogate the signal and reflect off the target
            received_target_sig = obj.channel_target(received_target_sig,defender_pos,tgt_pos,defender_vel,tgt_vel);
            received_target_sig = obj.SimulatedTarget.radar_target(received_target_sig);
            
            %receive the signal
            received_target_sig = obj.Defender.receiver(received_target_sig);
        end

        function received_attacker_sig = FMCW_simulate_received_chirp_from_attacker (obj,defender_pos, defender_vel, attacker_pos,attacker_vel,attacker_sig)
            %{
                Purpose: Simulate the signal that the defender would
                    receive from an interfering attacker
                Inputs:
                    [defender_pos, defender_vel]: the position and velocity
                        vectors for the defender (generated from 
                        phased.platform)
                    [attacker_pos,attacker_vel]: the position and velocity vectors
                        for the attacker (generated from phased.platform)
                    attacker_sig: the sig values of the
                        attacker's FMCW chirp waveform
                Outputs:
                    received_attacker_sig: the FMCW waveform signal after it has been
                        received by the defender's receiver, but before
                        its been dechirped
            %}
            
            %transmit the defender chirp
            received_attacker_sig = obj.Attacker.transmitter(attacker_sig);

            %propogate the signal and reflect off the target
            received_attacker_sig = obj.channel_attacker(received_attacker_sig,attacker_pos,defender_pos,attacker_vel,defender_vel);
            
            %receive the signal
            received_attacker_sig = obj.Defender.receiver(received_attacker_sig);
        end
        
        function sampled_IF_sig = FMCW_simulate_mixer_and_adc_sampling(obj,defender_sig,received_sig)
            %{
                Purpose: simulates a received signal going through a mixer
                    and then getting sampled by the defender's ADC at its ADC
                    sampling frequency. Returns the sampled IF that would
                    normally be recorded by the DCA1000
                Inputs:
                    defender_sig: the waveform of the signal transmitted by
                        the defender
                    received_sig: the signal that's received by the
                        defender after reflecting off of any targets. Includes
                        any attacker interference as well
                Outputs:
                    sampled_IF_sig: the sampled IF signal that would be
                        recorded by the DCA1000. Represents all of the samples
                        recorded for the given defender chirp
            %}
            %dechirp the received signal
            sampled_IF_sig = dechirp(received_sig,defender_sig);

            %put the sampled_IF_signal through the lowpass filter to simulate the mixer
            %and ADC
            sampled_IF_sig = filter(obj.Defender.low_pass_filter,sampled_IF_sig);
            
            %select only the samples that would have fallen in the sampling period
            num_samples_prior_to_sample_period = int32((obj.Defender.ADC_ValidStartTime_us * 1e-6) * obj.Defender.FMCW_sampling_rate_Hz);
            num_samples_in_sampling_period = obj.Defender.ADC_Samples * obj.Defender.downsample_factor;
            sampled_IF_sig = sampled_IF_sig(num_samples_prior_to_sample_period:num_samples_prior_to_sample_period + num_samples_in_sampling_period -1);
            
            %downsample the received signal to the sampling frequency of the receiver
            sampled_IF_sig = downsample(sampled_IF_sig,obj.Defender.downsample_factor);
        end
    
%% [3.3] FMCW simulations - compute IF signal for a given chirp

        function [sampled_IF_sig,combined_received_sigs, defender_waveform_sig, defender_waveform_t] = FMCW_simulate_scenario_chirp(obj,defender_chirp,defender_frame)
            %{
                Purpose: This function simulates a full FMCW scenario for a
                    specified chirp. For the specified chirp, we generate
                    the defender waveform, identify and generate the
                    relevant attacker chirps, simulate a target, transmit
                    the defender's chirp, reflect it off of the target,
                    determine what is received by the defender, and finally
                    output the IF signal which is identical to what is
                    generated by the actual hardware
                Inputs:
                    defender_chirp: the defender chirp to be simulated
                    defender_frame: the defender frame to be simulated
                Outputs:
                    sampled_IF_sig: the IF signal generated by the
                        defender. It is identical to the complex IF signal
                        measured by a physical radar in terms of number of
                        samples and sampling rate
                    combined_receive_sigs: this is a FMCW waveform signal
                        that represents the combination of any signal
                        received by the defender radar. It includes
                        reception from a target and any interfering
                        attacker chirps.
                    [defender_waveform_sig,defender_waveform_t]: the
                        [t,sig] values corresponding to the defender's
                        generated chirp
            %}

            %compute the target, defender, and attacker positions and
            %velocities
            [defender_pos, defender_vel,attacker_pos, attacker_vel, tgt_pos,tgt_vel] = obj.FMCW_determine_positions_and_velocities(defender_frame,defender_chirp);

            %generate the defender waveform
            [defender_waveform_t,defender_waveform_sig, defender_chirp_start_us,defender_chirp_end_us] = obj.FMCW_generate_defender_chirp_waveform(defender_frame,defender_chirp);
            
            %simulate the received signal from the target
            target_sig = obj.FMCW_simulate_received_chirp_from_target(defender_pos,defender_vel,tgt_pos,tgt_vel, defender_waveform_sig);
            
            %simulate the received signal from the attacker
            [attacker_t,attacker_sig] = obj.FMCW_generate_attacker_chirp_waveforms(defender_waveform_t,defender_waveform_sig,defender_chirp_start_us,defender_chirp_end_us);
            received_attacker_sig = obj.FMCW_simulate_received_chirp_from_attacker(defender_pos, defender_vel, attacker_pos,attacker_vel,attacker_sig);
            
            %combine the received signals from the target and the
            %interfering attacker chirps
            if obj.Attacker.enabled == true && obj.SimulatedTarget.enabled == true
                combined_received_sigs = target_sig + received_attacker_sig;
            elseif obj.Attacker.enabled == true
                combined_received_sigs = received_attacker_sig;
            elseif obj.SimulatedTarget.enabled == true
                combined_received_sigs = target_sig;
            else
                combined_received_sigs = zeros(size(defender_waveform_sig,1),1);
            end
            %generate the sampled_IF_signal
            sampled_IF_sig = obj.FMCW_simulate_mixer_and_adc_sampling(defender_waveform_sig,combined_received_sigs);
        end
 
%% [3.4] FMCW simulations - generate plots based on the IF data

        function FMCW_plot_sampled_IF_data(obj,axis, radar, sampled_IF_sig)
            %{
                Purpose: plots the real and imaginary samples for a given
                    sampled IF signal
                Inputs:
                    axis: the axis to plot to
                    radar: the radar object that has captured the IF signal
                    sampled_IF_sig: the sampled signal to be plotted
            %}
            sampling_period = 1/radar.ADC_SampleRate_MSps;
            t = 0:sampling_period:sampling_period*radar.ADC_Samples - sampling_period;
            real_data = real(sampled_IF_sig);
            imaginary_data = imag(sampled_IF_sig);
            plot(axis,t,real_data,t,imaginary_data);
            xlabel('time (us)');
            ylabel('codes');
            title_str = sprintf('ADC data');
            title(title_str);
            legend('real','imaginary');
        end
        
        function [fft_out,ranges] = FMCW_compute_range_fft(obj,radar,sampled_IF_sig)
            %{
                Purpose: computes the fft for a given sampled IF signal
                Inputs:
                    radar: the radar object that has captured the IF signal
                    sampled_IF_sig: the sampled signal to be plotted
                Outputs:
                    fft_out: the computed fft values (magnitude only)
                    ranges: all of the range bins for the computed fft
            %}
            fft_out = 20 * log10(abs(fft(sampled_IF_sig))); 
            ranges = (0:radar.ADC_Samples - 1) * radar.Range_Res_m;
        end
        
        function FMCW_plot_range_fft(obj,axis,radar,sampled_IF_sig)
            %{
                Purpose: plots the fft for a given sampled IF signal
                Inputs:
                    axis: the axis to plot to
                    radar: the radar object that has captured the IF signal
                    sampled_IF_sig: the sampled signal to be plotted
            %}
            [fft_out,ranges] = obj.FMCW_compute_range_fft(radar,sampled_IF_sig);
            plot(axis,ranges,fft_out);
            xlabel('range (m)');
            ylabel('magnitude (dB)');
            title_str = sprintf('Range FFT ');
            title(title_str);
        end
    end
    methods (Access = private)

        function plot_radar_chirp(obj,axis,radar, tx,f_tx,t_sampling,f_sampling)
            %{
            Purpose: For a given radar, plot a chirp starting at a
                specified time
            Inputs:
                axis: the axis variable for the figure to plot in
                radar: the radar to plot from
                chirp_start_time_us: the start time of the chirp to be
                    plotted
                [tx,f_tx]: the t,f values to be plotted for the chirp
                [t_sample,f_sampling]: the t,f values for the sampling
                    period to be plotted for the chirp
            %}
            

            %plot f over t for the whole Tx On period
            p = scatter(axis,tx,f_tx,'.');
            p.DisplayName = sprintf('%s Tx On',radar.radar_name);
            p.MarkerEdgeColor = radar.tx_period_plot_color;
            
            %plot f over t the Tx sampling period
            hold(axis,'on');
            
            p = scatter(axis,t_sampling,f_sampling,'.');
            p.DisplayName = sprintf('%s Sampling Period',radar.radar_name);
            p.MarkerEdgeColor = radar.tx_sampling_period_plot_color;
            
            title('Chirp frequency over time')
            xlabel('time (us)')
            ylabel('frequency (GHz)')
            legend;
            
            hold(axis,'off');
        end
    
        function plot_radar_frame(obj,axis,radar, tx,f_tx,t_sampling,f_sampling)
            %{
                Purpose: For a given radar, plot a chirp starting at a
                    specified time
                Inputs:
                    axis: the axis variable for the figure to plot in
                    radar: the radar to plot from
                    frame_start_time_us: the start time of the frame to be
                        plotted
                    [tx,f_tx]: the t,f values to be plotted for the frame
                    [t_sample,f_sampling]: the t,f values for the sampling
                        period to be plotted for the frame
            %}
            
            %plot f over t for the whole Tx On period
            p = scatter(axis,tx,f_tx,'.');
            p.DisplayName = sprintf('%s Tx On',radar.radar_name);
            p.MarkerEdgeColor = radar.tx_period_plot_color;
            
            %plot f over t the Tx sampling period
            hold(axis,'on')
            p = scatter(axis,t_sampling,f_sampling,'.');
            p.DisplayName = sprintf('%s Sampling Period',radar.radar_name);
            p.MarkerEdgeColor = radar.tx_sampling_period_plot_color;
            
            title('Frame frequency over time')
            xlabel('time (ms)')
            ylabel('frequency (GHz)')
            legend;
            hold(axis,'off');
        end
    
        function [attack_region_t,attack_region_f] = identify_valid_chirp_interference(obj, defender_t_sampling, defender_f_sampling, attacker_tx, attacker_f_tx)
            %{
                Purpose: identifies the points in the attacker chirp that
                    would actually interfere with the defender chirp (i.e: the
                    difference between the attacker and defender chirps at a
                    specific value of t is within half of the sampling
                    frequency of the defender
                Inputs:
                    [defender_t_sampling,defender_f_sampling]: the [t,f]
                        values of the sampling period for the defender
                    [attacker_tx, attacker_f_tx]: the [t,f] values of the
                        sampling period for the transmission period of the
                        attacker
                Outputs:
                    [attack_region_t, attack_region_f]: the [t,f] values
                    corresponding to the when the attacker is actually
                    interfering with the defender
            %}

            %identify intersection between defender_t_sampling and attacker_tx
            %appears to be a floating point discrepancy in the time values
            %and so try and reset them so that they are equal again

            defender_t_sampling = round(defender_t_sampling/obj.Defender.plotResolution_us) * obj.Defender.plotResolution_us;
            attacker_tx = round(attacker_tx / obj.Defender.plotResolution_us) * obj.Defender.plotResolution_us;

            [intersection_values,i_defender_t_sampling,i_attacker_tx] = intersect(defender_t_sampling,attacker_tx);
            intersection_defender_t = defender_t_sampling(i_defender_t_sampling);
            intersection_defender_f = defender_f_sampling(i_defender_t_sampling);
    
            intersection_attacker_t = attacker_tx(i_attacker_tx);
            intersection_attacker_f = attacker_f_tx(i_attacker_tx);
            
            %determine the maximum if frequency observable from the defender
            max_if_freq_GHz = (obj.Defender.ADC_SampleRate_MSps/2)/1000;
            valid_intersection = (abs(intersection_attacker_f - intersection_defender_f) < max_if_freq_GHz);
    
            attack_region_t = intersection_attacker_t(valid_intersection);
            attack_region_f = intersection_attacker_f(valid_intersection);
        end

        function [attack_region_t,attack_region_f] = identify_valid_frame_interference(obj, defender_t_sampling, defender_f_sampling, attacker_tx, attacker_f_tx)
            %{
                Purpose: identifies the points in the attacker chirp that
                    would actually interfere with the defender chirp (i.e: the
                    difference between the attacker and defender chirps at a
                    specific value of t is within half of the sampling
                    frequency of the defender
                Inputs:
                    [defender_t_sampling,defender_f_sampling]: the [t,f]
                        values of the sampling period for the defender
                    [attacker_tx, attacker_f_tx]: the [t,f] values of the
                        sampling period for the transmission period of the
                        attacker
                Outputs:
                    [attack_region_t, attack_region_f]: the [t,f] values
                    corresponding to the when the attacker is actually
                    interfering with the defender
            %}

            %identify intersection between defender_t_sampling and attacker_tx
            %appears to be a floating point discrepancy in the time values
            %and so try and reset them so that they are equal again
            
            %multiplying and dividing by 1000 to convert the resolution to
            %ms
            
            defender_t_sampling = round(defender_t_sampling*1000/obj.Defender.plotResolution_us) * obj.Defender.plotResolution_us/1000;
            attacker_tx = round(attacker_tx *1000 / obj.Defender.plotResolution_us) * obj.Defender.plotResolution_us/1000;

            [intersection_values,i_defender_t_sampling,i_attacker_tx] = intersect(defender_t_sampling,attacker_tx);
            intersection_defender_t = defender_t_sampling(i_defender_t_sampling);
            intersection_defender_f = defender_f_sampling(i_defender_t_sampling);
    
            intersection_attacker_t = attacker_tx(i_attacker_tx);
            intersection_attacker_f = attacker_f_tx(i_attacker_tx);
            
            %determine the maximum if frequency observable from the defender
            max_if_freq_GHz = (obj.Defender.ADC_SampleRate_MSps/2)/1000;
            valid_intersection = (abs(intersection_attacker_f - intersection_defender_f) < max_if_freq_GHz);
    
            attack_region_t = intersection_attacker_t(valid_intersection);
            attack_region_f = intersection_attacker_f(valid_intersection);
        end

        function [intersection_attacker_t,intersection_attacker_sig] = FMCW_identify_valid_chirp_intersection(obj, defender_t, defender_sig, attacker_t, attacker_sig)
            %{
                Purpose: identifies points in the attacker chirp that
                    would fall within the transmit period of the defender
                Inputs:
                    [defender_t,defender_sig]: the [t,sig]
                        values of the defender
                    [attacker_t, attacker_sig]: the [t,sig] values of
                    attacker chirp
                Outputs:
                    [intersection_attacker_t,intersection_attacker_sig]: the [t,sig]
                    values corresponding to when the attacker chirp falls
                    within the defender chirp
                Note: Times are in us, not seconds
            %}

            
            %identify the start and stop index (in the attacker array) of
            %the overlap (if there is any overlap at all)
            tolerance = 1/obj.Defender.FMCW_sampling_rate_Hz * 1e6;
            attacker_start_index = 0;
            attacker_end_index = 0;

            defender_start_index = 0;
            defender_end_index = 0;

            %solve for the start index first
            [M,I] = min(abs(attacker_t(1) - defender_t));
            if M < tolerance %see if its the first value of the attacker
                attacker_start_index = 1;
                defender_start_index = I;
            else %if first time value isn't in the attacker chirp, check the first value of the defender
                [M,I] = min(abs(defender_t(1) - attacker_t));
                if M < tolerance
                    attacker_start_index = I;
                    defender_start_index = 1;
                else %at this point, the attacker chirp must not overlap the defender chirp
                    attacker_start_index = 0;
                    defender_start_index = 0;
                end
            end

            %solve for the end index next
            if attacker_start_index >= 1
                [M,I] = min(abs(attacker_t(end) - defender_t));
                if M < tolerance %check if its the last value of the attacker
                    attacker_end_index = size(attacker_t,1);
                    defender_end_index = I;
                else % must be the the last value of the defender
                    [M,I] = min(abs(defender_t(end) - attacker_t));
                    attacker_end_index = I;
                    defender_end_index = size(defender_t,1);
                end
            else
                attacker_end_index = 0;
                defender_end_index = 0;
            end

            %once the overlap has been identified, filter the attacker
            %arrays to only be the part that has overlap
            
            if attacker_start_index >= 1 && attacker_end_index >= 1
                intersection_attacker_t = defender_t(defender_start_index:defender_end_index);
                intersection_attacker_sig = attacker_sig(attacker_start_index:attacker_end_index);
            else
                intersection_attacker_t = [];
                intersection_attacker_sig = [];
            end
        end
    end
end