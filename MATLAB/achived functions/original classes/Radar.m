%% ARCHIVED Verion DO NOT USE THIS CODE

classdef Radar < handle
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here

    properties (Access = public)
        %chirp parameters
        StartFrequency_GHz
        FrequencySlope_MHz_us
        IdleTime_us
        TxStartTime_us
        ADC_ValidStartTime_us
        ADC_Samples
        ADC_SampleRate_MSps
        RampEndTime_us
        Chirp_Tx_Bandwidth_MHz          %calculated value, not setable in UI
        Chirp_Sampling_Bandwidth_MHz    %calculated value, not setable in UI
        ADC_SamlingPeriod_us            %calculated value, not setable in UI
        ChirpCycleTime_us               %calculated value, not setable in UI
        Lambda_m                        %calculated value, not setable in UI
        
        %frame parameters
        NumChirps
        FramePeriodicity_ms
        ActiveFrameTime_ms              %calculated value, not setable in UI

        %radar performance specs        %calculated values, not setable in UI
        Range_Max_m
        Range_Res_m
        V_Max_m_per_s
        V_Res_m_per_s

        %parameters for the position and speed of the radar
        position_m
        velocity_m_per_s
        
        %if the radar is an attacker, use this parameter to set its offset
        %in us from the defender
        timing_offset_us

        %parameter to enable or disable a given radar (use for enabling or
        %disabling the attacker)
        enabled
        
        %% The following parameters don't appear in the UI
        %define parameters for plotting
        plotResolution_us
        tx_period_plot_color
        tx_sampling_period_plot_color
        radar_name

        %parameters for the FMCW waveform, not setable parameters in UI
        waveform
        FMCW_sampling_rate_Hz
        downsample_factor
        sweep_time
        additional_BW_MHz                   %for if the start frequency isn't 77.0 GHz
        additional_sweep_time_us            %for if the start frequency isn't 77.0 GHz

        %parameters for the transmitter and receiver
        transmitter         %phased.Transmitter object
        receiver            %phased.Receiver

        %parameter for the platform of the radar
        platform            %phased.Platform object

        %parameter for a filter to simulate ADC attenuation of higher
        %frequency elements
        low_pass_filter

    end
    properties(Access = private)
        %parameters for the Rx and Tx power when dealing with the FMCW
        %waveform, not setable parameters in the UI
        ant_aperture_m2
        ant_gain_dB

        tx_power_W
        tx_gain_dB

        rx_gain_dB
        rx_nf_dB

        %parameters for the lowpass filter design
        Fp      %start frequency of pass band
        Fst     %start ofstop pand 
        Ap      %ripple to allow in the pass band
        Ast     %attenuation in the stop band
        Fs      %sampling frequency
    end

    methods (Access = public)
        function obj = Radar()
            %UNTITLED3 Construct an instance of this class
            %   Detailed explanation goes here
        end

        function print_chirp_parameters(obj)
            fprintf("Chirp Parameters\n")
            fprintf("\t Start Frequency: \t\t %.2f GHz\n",obj.StartFrequency_GHz)
            fprintf("\t Frequency Slope: \t\t %.2f MHz/us\n",obj.FrequencySlope_MHz_us)
            fprintf("\t Idle Time: \t\t\t %.2f us\n",obj.IdleTime_us)
            fprintf("\t Tx Start Time: \t\t %.2f us\n",obj.TxStartTime_us)
            fprintf("\t ADC Valid Start Time: \t\t %.2f us\n",obj.ADC_ValidStartTime_us)
            fprintf("\t ADC Samples: \t\t\t %d \n",obj.ADC_Samples)
            fprintf("\t ADC Sample Rate: \t\t %.2f MSps\n",obj.ADC_SampleRate_MSps)
            fprintf("\t Ramp End Time: \t\t %.2f us\n",obj.RampEndTime_us)
            fprintf("\t Chirp Tx Bandwidth: \t\t %.2f MHz\n",obj.Chirp_Tx_Bandwidth_MHz)
            fprintf("\t Chirp Sampling Bandwidth: \t %.2f MHz\n", obj.Chirp_Sampling_Bandwidth_MHz)
            fprintf("\t ADC Sampling Period: \t\t %.2f us\n",obj.ADC_SamlingPeriod_us)
            fprintf("\t Chirp Cycle Time: \t\t %.2f us\n",obj.ChirpCycleTime_us)
            fprintf("\t Chirp Wavelength: \t\t %.2f nm\n",obj.Lambda_m * 1e-9)
        end

        function print_frame_parameters(obj)
            fprintf("Frame Parameters\n")
            fprintf("\t Number of Chirps \t\t %d\n",obj.NumChirps)
            fprintf("\t FramePeriodicity \t\t %0.2f ms\n",obj.FramePeriodicity_ms)
            fprintf("\t Active Frame Time \t\t %0.2f ms\n",obj.ActiveFrameTime_ms)
        end
        
        function print_performance_specs(obj)
            fprintf("Performance Specifications\n")
            fprintf("\t Max Range \t\t\t %0.2f m\n",obj.Range_Max_m)
            fprintf("\t Range Resolution \t\t %0.2f m\n",obj.Range_Res_m)
            fprintf("\t Max Velocity \t\t\t %0.2f m/s\n",obj.V_Max_m_per_s)
            fprintf("\t Velocity Resolution \t\t %0.2f m/s\n",obj.V_Res_m_per_s)
        end

        function [t,f] = generate_chirp_f_over_t_vals(obj,chirp_start_time_us)
%{
            Purpose: generates (t,f) coordinates for the Tx
                    frequency over time
            Inputs:
                    chirp_start_time: time that a chirp starts at
            Outputs:
                   t,f: (t,f) coordinates for the Tx frequency over time
%}            

            %round the times to be within the level of precision
            chirp_cycle_time = round((obj.ChirpCycleTime_us + chirp_start_time_us)/ obj.plotResolution_us)*obj.plotResolution_us;
            idle_time = round((obj.IdleTime_us + chirp_start_time_us) / obj.plotResolution_us)*obj.plotResolution_us;
            tx_start_time = round(obj.TxStartTime_us/ obj.plotResolution_us)*obj.plotResolution_us;
            
            %calculate the sloping values
            t = idle_time:obj.plotResolution_us:chirp_cycle_time;
            f = obj.StartFrequency_GHz + obj.FrequencySlope_MHz_us * (t - obj.IdleTime_us - chirp_start_time_us)/1000;

            %If Tx Start Time is not zero
            if obj.TxStartTime_us ~= 0
                t_additional = idle_time - tx_start_time:obj.plotResolution_us:idle_time - obj.plotResolution_us;
                t = [t_additional, t];

                f_additional =  obj.StartFrequency_GHz * ones(1,size(t_additional,2));
                f = [f_additional, f];
            end
        end

        function [ts,f] = generate_chirp_f_over_Tsampling_vals(obj, chirp_start_time_us)
%{
            Purpose: generates (t,f) coordinates for the Tx frequency over
                time for a single chirp
            Inputs:
                    chirp_start_time: time that a chirp starts at
            Outputs:
                t,f: (t,f) coordinates for the Tx frequency over the
                    sampling period
%}
            adc_sampling_start_time = round((obj.IdleTime_us + obj.ADC_ValidStartTime_us + chirp_start_time_us)/obj.plotResolution_us)*obj.plotResolution_us;
            adc_sampling_period = round(obj.ADC_SamlingPeriod_us/obj.plotResolution_us) * obj.plotResolution_us;

            %calculate the currently sloping values
            ts = adc_sampling_start_time:obj.plotResolution_us:adc_sampling_start_time + adc_sampling_period;
            f = obj.StartFrequency_GHz + obj.FrequencySlope_MHz_us * (ts - obj.IdleTime_us - chirp_start_time_us)/1000;
        end

        function [t,f] = generate_frame_f_over_T_vals(obj,frame_start_time_ms)
            %{
            Purpose: generates (t,f) coordinates for the Tx frequency over
                    time for an entire frame (in ms)
            Outputs:
                t,f: (t,f) coordinates for the Tx frequency over the
                    period of each chirp in the frame (in ms)
            %}
            
            %pre-allocate the size of the output array

            %round the times to be within the level of precision
            chirp_cycle_time = round((obj.ChirpCycleTime_us)/ obj.plotResolution_us)*obj.plotResolution_us;
            idle_time = round((obj.IdleTime_us) / obj.plotResolution_us)*obj.plotResolution_us;
            tx_start_time = round(obj.TxStartTime_us/ obj.plotResolution_us)*obj.plotResolution_us;
            
            %calculate the sloping values
            t_size = int32((chirp_cycle_time - idle_time)/obj.plotResolution_us + 1);

            %If Tx Start Time is not zero
            if obj.TxStartTime_us ~= 0
                t_size = t_size + tx_start_time/obj.plotResolution_us;

            end
            
            t = zeros(1,t_size * obj.NumChirps);
            f = zeros(1,t_size * obj.NumChirps);
            next_start_index = 1;

            %generate the f and t arrays 
            for chirp = 1:obj.NumChirps
                chirp_start_time_us = obj.ChirpCycleTime_us * (chirp - 1) + (frame_start_time_ms * 1000);
                [t_chirp,f_chirp] = obj.generate_chirp_f_over_t_vals(chirp_start_time_us);

                if(size(t_chirp,2) ~= t_size)
                    difference = size(t_chirp,2) - t_size;
                    if(difference > 0)
                        t = [t,zeros(1,difference)];
                        f = [f,zeros(1,difference)];
                    else
                        t = t(1:end + difference); %will reduce its size
                        f = f(1:end + difference); %will reduce its size
                    end
                    
                    t(next_start_index: next_start_index + t_size + difference - 1) = t_chirp;
                    f(next_start_index: next_start_index + t_size + difference - 1) = f_chirp;
                    next_start_index = next_start_index + t_size + difference;
                else
                    t(next_start_index: next_start_index + t_size - 1) = t_chirp;
                    f(next_start_index: next_start_index + t_size - 1) = f_chirp;
                    next_start_index = next_start_index + t_size;
                end
            end
            t = t/1000;
        end

        function [t,f] = generate_frame_f_over_Tsampling_vals(obj, frame_start_time_ms)
            %{
            Purpose: generates (t,f) coordinates for the Tx frequency over
                    the sampling period for each chirp in an entire frame
                    (in ms)
            Outputs:
                t,f: (t,f) coordinates for the Tx frequency over the
                    sampling period of each chirp in the frame (in ms)
            %}
            
            %pre-allocate the size of the output array

            %round the times to be within the level of precision
            adc_sampling_period = round(obj.ADC_SamlingPeriod_us/obj.plotResolution_us) * obj.plotResolution_us;
            
            %calculate the sloping values
            t_size = int32((adc_sampling_period)/obj.plotResolution_us + 1);
            
            t = zeros(1,t_size * obj.NumChirps);
            f = zeros(1,t_size * obj.NumChirps);
            next_start_index = 1;

            %generate the f and t arrays 
            for chirp = 1:obj.NumChirps
                chirp_start_time_us = obj.ChirpCycleTime_us * (chirp - 1) + (frame_start_time_ms * 1000);
                [t_chirp,f_chirp] = obj.generate_chirp_f_over_Tsampling_vals(chirp_start_time_us);

                if(size(t_chirp,2) ~= t_size)
                    difference = size(t_chirp,2) - t_size;
                    if(difference > 0)
                        t = [t,zeros(1,difference)];
                        f = [f,zeros(1,difference)];
                    else
                        t = t(1:end + difference); %will reduce its size
                        f = f(1:end + difference);
                    end
                    
                    t(next_start_index: next_start_index + t_size + difference - 1) = t_chirp;
                    f(next_start_index: next_start_index + t_size + difference - 1) = f_chirp;
                    next_start_index = next_start_index + t_size + difference;
                else
                    t(next_start_index: next_start_index + t_size - 1) = t_chirp;
                    f(next_start_index: next_start_index + t_size - 1) = f_chirp;
                    next_start_index = next_start_index + t_size;
                end
            end
            t = t/1000;
        end

        function compute_calculated_vals(obj)
            %{
                Purpose: computes all of the Radar class parameters that
                are marked as "calculated value, not setable in UI"
            %}

            %compute calculated chirp parameters
            obj.ADC_SamlingPeriod_us            = obj.ADC_Samples / obj.ADC_SampleRate_MSps;
            obj.Chirp_Tx_Bandwidth_MHz          = obj.FrequencySlope_MHz_us * obj.RampEndTime_us;
            obj.Chirp_Sampling_Bandwidth_MHz    = obj.FrequencySlope_MHz_us * obj.ADC_SamlingPeriod_us;
            obj.ChirpCycleTime_us               = obj.IdleTime_us + obj.RampEndTime_us;
            obj.Lambda_m                        = physconst('LightSpeed') / obj.StartFrequency_GHz * 1e-9;

            %compute calculated frame parameters
            obj.ActiveFrameTime_ms              = obj.ChirpCycleTime_us * obj.NumChirps/1000;

            %calculate radar performance parameters
            obj.Range_Max_m = physconst('LightSpeed') * obj.ADC_SampleRate_MSps*1e6 / (2*obj.FrequencySlope_MHz_us*1e6*1e6);
            obj.Range_Res_m = physconst('LightSpeed') / (2*obj.Chirp_Sampling_Bandwidth_MHz*1e6);
            obj.V_Max_m_per_s = obj.Lambda_m / (4* obj.ChirpCycleTime_us * 1e-6);
            obj.V_Res_m_per_s = obj.Lambda_m / (2*obj.NumChirps*obj.ChirpCycleTime_us * 1e-6);
        end

        function configure_transmitter_and_receiver(obj)
            %{
                Purpose: configures the transmitter and receiver objects
                Note: all other parameters must be set before configuring
                the transmitter and receiver
            %}
            lambda = freq2wavelen(obj.StartFrequency_GHz * 1e9);


            obj.ant_aperture_m2 = 6.06e-4;                                  % in square meter
            obj.ant_gain_dB = aperture2gain(obj.ant_aperture_m2,lambda);    % in dB
            
            obj.tx_power_W = db2pow(5)*1e-3;                                % in watts
            obj.tx_gain_dB = 9+ obj.ant_gain_dB;                            % in dB
            
            obj.rx_gain_dB = 15+ obj.ant_gain_dB;                           % in dB
            obj.rx_nf_dB = 4.5;                                             % in dB
            
            obj.transmitter = phased.Transmitter( ...
                'PeakPower',obj.tx_power_W, ...
                'Gain',obj.tx_gain_dB);
            obj.receiver = phased.ReceiverPreamp( ...
                'Gain',obj.rx_gain_dB, ...
                'NoiseFigure',obj.rx_nf_dB,...
                'SampleRate',obj.FMCW_sampling_rate_Hz);
        end

        function configure_lowpass_filter(obj)
            %{
                Purpose: function configures the lowpass filter that will
                be used to simulate how the ADC and mixer attenuate higher
                frequency components in a signal
            %}
            
            obj.Fp = 20 * 1e6;      %start frequency of pass band
            obj.Fst = 30 * 1e6;     %start of stop band 
            obj.Ap = 0.5;           %ripple to allow in the pass band
            obj.Ast = 40;           %attenuation in the stop band
            obj.Fs = obj.FMCW_sampling_rate_Hz; %sampling frequency
            
            d = fdesign.lowpass('Fp,Fst,Ap,Ast',obj.Fp,obj.Fst,obj.Ap,obj.Ast,obj.Fs);
            obj.low_pass_filter = design(d,'butter','MatchExactly','passband');
            %Hd = design(d,'equiripple');
            %fvtool(Hd)
        end

        function [t,sig] = FMCW_generate_chirp_sig_vals(obj,chirp_start_time_us)
            %{
                Purpose: computes the output signal of the radar's waveform
                    and the times that correspond to a specific chirp
                Inputs:
                        chirp_start_time: time that a chirp starts at in us
                Outputs:
                       t,sig: the output of from calling the radar's
                           waveform function and the times associated with a
                           specific chirp
            %}
            sig = obj.waveform();

            %if the starting frequency isn't 77Ghz, then we need to remove
            %the first part of the signal so that the signal starts at the
            %right frequency. This is done because the FMCWwaveform
            %function only allows us to start at 0Ghz, not a custom
            %frequency
            
            if obj.StartFrequency_GHz ~= 77.0
                samples_to_remove = int32(obj.additional_sweep_time_us * 1e-6 * obj.FMCW_sampling_rate_Hz);
                sig = sig(samples_to_remove + 1:end);
            end
            
            start_time = chirp_start_time_us + obj.IdleTime_us;
            increment = (1/obj.FMCW_sampling_rate_Hz) * 1e6;
            end_time = start_time + obj.sweep_time * 1e6 - obj.additional_sweep_time_us;
            t = start_time:increment:end_time - increment;
            t = reshape(t,size(t,2),1);
        end
    end
end