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
%}

classdef Simulator < handle
    %SIMULATOR class used to support simulation in the simulink model

    properties (Access = public)
        Attacker
        Victim
        SimulatedTarget

        %variables used in computing the FMCW simulation - MOVE TO PRIVATE
        channel_target          %phased.FreeSpace object for a target
        channel_attacker        %phased.FreeSpace object for an attacker
    end

    methods (Access = public)
        function obj = Simulator()
            %Simulator Construct an instance of this class
            %   Detailed explanation goes here
            obj.Attacker = Radar();
            obj.Victim = Radar();
        end

%% [1] Functions to configure the simulation
        
        function configure_FMCW_parameters(obj)
            %{
                Purpose: initializes all of the needed waveform parameters
                    to simulate the actual waveforms for the attacker and
                    victim
            %}

            %obj.Victim.downsample_factor = ceil((2 * obj.Victim.Chirp_Tx_Bandwidth_MHz * 1e6) / (obj.Victim.ADC_SampleRate_MSps * 1e6));
                %NOTE: if needed, can use nextpow2() function when
                %implementing on the fpga

            obj.Victim.downsample_factor = ceil((obj.Victim.Chirp_Tx_Bandwidth_MHz * 1e6) / (obj.Victim.ADC_SampleRate_MSps * 1e6));
            obj.Victim.decimation_factor = obj.Victim.downsample_factor;
            
            %removed the doubling of the chirp Tx Bandwidth since we are
            %assuming complex sampling

            obj.Victim.FMCW_sampling_rate_Hz = obj.Victim.ADC_SampleRate_MSps * 1e6 * obj.Victim.downsample_factor;
            obj.Attacker.FMCW_sampling_rate_Hz = obj.Victim.FMCW_sampling_rate_Hz;

            %set sweep time
            obj.Victim.sweep_time = obj.Victim.RampEndTime_us * 1e-6;
            obj.Attacker.sweep_time = obj.Attacker.RampEndTime_us * 1e-6;

            %configure the number of samples in the idle time and in each
            %frame
            
            %configure the FMCW waveforms for the victim and attacker
            %configure waveforms
            obj.Victim.configure_waveform_and_chirp();
            obj.Attacker.configure_waveform_and_chirp();
            
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

            obj.Attacker.configure_transmitter_and_receiver();
            %obj.Attacker.configure_lowpass_filter();
            %obj.Attacker.configure_CFAR_detector();

            obj.Victim.configure_transmitter_and_receiver();
            obj.Victim.configure_lowpass_filter();
            obj.Victim.configure_decimator();
            obj.Victim.configure_RangeDopplerResponse();
            obj.Victim.configure_CFAR_detector();
            obj.Victim.configure_DB_scan();

            %
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
            obj.Victim.NumChirps                  = 32;
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
            obj.Attacker.StartFrequency_GHz         = 77.0;
            obj.Attacker.FrequencySlope_MHz_us      = 10.76;
            obj.Attacker.TxStartTime_us             = 0;
            obj.Attacker.ADC_Samples                = 256;
            obj.Attacker.ADC_SampleRate_MSps        = 7.17;
            obj.Attacker.ChirpCycleTime_us          = 50;
            
            %setup the victim's frame parameters
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
        end

        function load_B210_victim_params(obj)
            %setup the victim's chirp parameters
            obj.Victim.StartFrequency_GHz         = 77.0;
            obj.Victim.FrequencySlope_MHz_us      = 1.2;
            obj.Victim.TxStartTime_us             = 0;
            obj.Victim.ADC_Samples                = 64;
            obj.Victim.ADC_SampleRate_MSps        = 1.54;
            obj.Victim.ChirpCycleTime_us          = 50;             

            %setup the victim's frame parameters
            obj.Victim.NumChirps                  = 32;
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
            obj.Attacker.StartFrequency_GHz         = 77.0;
            obj.Attacker.FrequencySlope_MHz_us      = 1.2;
            obj.Attacker.TxStartTime_us             = 0;
            obj.Attacker.ADC_Samples                = 64;
            obj.Attacker.ADC_SampleRate_MSps        = 1.54;
            obj.Attacker.ChirpCycleTime_us          = 50;             

            %setup the victim's frame parameters
            obj.Attacker.NumChirps                  = 32;
            obj.Attacker.FramePeriodicity_ms        = 33.33;
            
            %define plot color default values
            obj.Attacker.plotResolution_us = 0.01;
            obj.Attacker.tx_period_plot_color = 'red';
            obj.Attacker.tx_sampling_period_plot_color = 'magenta';
            obj.Attacker.radar_name = 'Attacker';

            %set timing offset to zero as this is the victim
            obj.Attacker.timing_offset_us = 0;

            %compute all remaining "calculated" values
            obj.Attacker.compute_calculated_vals();
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

            obj.SimulatedTarget = Target(position_m,velocity_meters_per_s,rcs_sq_meters,operating_frequency_Hz);
        end
        
        function load_realistic_attacker_and_victim_position_and_velocity(obj)
            %{
                Purpose: configures a default scenario for the attacker and
                victim positions and velocities
            %}
            obj.Victim.position_m = [0;0;0];
            obj.Victim.velocity_m_per_s = [0;0;0];
            obj.Victim.platform = phased.Platform( ...
                'InitialPosition',obj.Victim.position_m, ...
                'Velocity',obj.Victim.velocity_m_per_s);

            obj.Attacker.position_m = [75;0;0];
            obj.Attacker.velocity_m_per_s = [0;0;0];
            obj.Attacker.platform = phased.Platform( ...
                'InitialPosition',obj.Attacker.position_m, ...
                'Velocity',obj.Attacker.velocity_m_per_s);
        end
    
%% [2] Functions for running the FMCW Simulation on Matlab

function [victim_pos, victim_vel,attacker_pos, attacker_vel, tgt_pos,tgt_vel] = FMCW_determine_positions_and_velocities(obj,victim_frame,victim_chirp)
        %{
            Purpose: determine the position and velocity for the
                attacker,defender, and target at the start of a given chirp
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
            %calculate the time that the desired chirp will start at
            frame_start_time_s = (obj.Victim.FramePeriodicity_ms * 1e-3) * (victim_frame - 1);
            chirp_start_time_s = frame_start_time_s + (obj.Victim.ChirpCycleTime_us * 1e-6) * (victim_chirp -1);
            
            %reset the positions of each platform to be safe
            reset(obj.Victim.platform);
            reset(obj.Attacker.platform);
            reset(obj.SimulatedTarget.platform);

            if chirp_start_time_s <= 0
                time_increment = obj.Victim.ChirpCycleTime_us * 1e-6;
                [victim_pos, victim_vel] = obj.Victim.platform(time_increment);
                [tgt_pos,tgt_vel] = obj.SimulatedTarget.platform(time_increment);
                [attacker_pos, attacker_vel] = obj.Attacker.platform(time_increment);
            else
                %take a step for each platform so that they update to the
                %values at the start of the frame
                obj.Victim.platform(chirp_start_time_s);
                obj.SimulatedTarget.platform(chirp_start_time_s);
                obj.Attacker.platform(chirp_start_time_s);
    
                %repeat the process to obtain the positions and velocities at
                %the start of the given chirp
                [victim_pos, victim_vel] = obj.Victim.platform(chirp_start_time_s);
                [tgt_pos,tgt_vel] = obj.SimulatedTarget.platform(chirp_start_time_s);
                [attacker_pos, attacker_vel] = obj.Attacker.platform(chirp_start_time_s);
            end
        end
        
%% [3] Functions to generate frequency over time for a given chirp or frame

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

        function plot_attacker_and_victim_chirps(obj,axis,victim_frame,victim_chirp)
            %{
                Purpose: generate a plot for a specific victim chirp in a
                    given frame and overlay any attacker chirps that would
                    overlap with that victim chirp
                Inputs:
                    axis: the axis object corresponding to the figure to
                        plot in
                    victim_frame: the desired frame number of the
                        victim
                    victim_chirp: the desired chirp number
            %}

            %calculate the [t,f] values for the victim chirp

            victim_frame_start_us = obj.Victim.FramePeriodicity_ms * 1000 * (victim_frame - 1) + obj.Victim.timing_offset_us;
            victim_chirp_start_us = victim_frame_start_us + obj.Victim.ChirpCycleTime_us * (victim_chirp - 1);
            victim_chirp_end_us = victim_frame_start_us + obj.Victim.ChirpCycleTime_us * (victim_chirp);
            
            [victim_tx,victim_f_tx] = obj.Victim.generate_chirp_f_over_t_vals(victim_chirp_start_us);
            [victim_t_sampling,victim_f_sampling] = obj.Victim.generate_chirp_f_over_Tsampling_vals(victim_chirp_start_us);

            %now, calculate the [t,f] values for the attacker
            
            %determine the frames that will need to be plotted
            attacker_frame = floor((victim_chirp_start_us - obj.Attacker.timing_offset_us)/(obj.Attacker.FramePeriodicity_ms*1000))+1;
            attacker_frame_start_us = obj.Attacker.FramePeriodicity_ms * 1000 * (attacker_frame - 1) + obj.Attacker.timing_offset_us;


            %initialize the attacker chirp to start with
            attacker_chirp = floor((victim_chirp_start_us - attacker_frame_start_us)/ obj.Attacker.ChirpCycleTime_us)+1;
            if attacker_chirp > obj.Attacker.NumChirps
                attacker_chirp = obj.Attacker.NumChirps;
            end

            %calculate the [t,f] values for the first chirp
            attacker_chirp_start_us = attacker_frame_start_us + obj.Attacker.ChirpCycleTime_us * (attacker_chirp - 1);
            attacker_chirp_end_us = attacker_frame_start_us + obj.Attacker.ChirpCycleTime_us * (attacker_chirp);
            [attacker_tx,attacker_f_tx] = obj.Attacker.generate_chirp_f_over_t_vals(attacker_chirp_start_us);
            [attacker_t_sampling,attacker_f_sampling] = obj.Attacker.generate_chirp_f_over_Tsampling_vals(attacker_chirp_start_us);

            while attacker_chirp_end_us < victim_chirp_end_us
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

            %identify the points where the victim would actually observe the
            
           [attack_region_t, attack_region_f] = obj.identify_valid_chirp_interference(victim_t_sampling, victim_f_sampling, attacker_tx, attacker_f_tx);

            %plot everything
            obj.plot_radar_chirp(axis,obj.Victim,victim_tx,victim_f_tx,victim_t_sampling,victim_f_sampling);
            hold on
            obj.plot_radar_chirp(axis,obj.Attacker,attacker_tx,attacker_f_tx,attacker_t_sampling,attacker_f_sampling);
            
            hold on
            p = scatter(axis,attack_region_t,attack_region_f,'.');
            p.DisplayName = sprintf('Valid attack region ');
            p.MarkerEdgeColor = 'green';
            hold off
            
            % by default, zoom the axis to only focus on the victim chirp
            set(axis,'XLim',[victim_chirp_start_us - 10,victim_chirp_end_us + 10]);
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

        function plot_attacker_and_victim_frames(obj,axis,victim_frame)
                %{
                Purpose: generate a plot for a specific victim frame and 
                    overlay any attacker frames that would overlap with 
                    that victim chirp
                Inputs:
                    axis: the axis object corresponding to the figure to
                        plot in
                    victim_frame: the desired frame number of the
                        victim
                %}
             %calculate the [t,f] values for the victim chirp

            victim_frame_start_ms = obj.Victim.FramePeriodicity_ms * (victim_frame - 1) + obj.Victim.timing_offset_us / 1000;
            victim_frame_end_ms = victim_frame_start_ms + obj.Victim.FramePeriodicity_ms * (victim_frame);
            
            [victim_tx,victim_f_tx] = obj.Victim.generate_frame_f_over_T_vals(victim_frame_start_ms);
            [victim_t_sampling,victim_f_sampling] = obj.Victim.generate_frame_f_over_Tsampling_vals(victim_frame_start_ms);

            %now, calculate the [t,f] values for the attacker
            
            %determine the frames that will need to be plotted
            attacker_frame = floor((victim_frame_start_ms - (obj.Attacker.timing_offset_us/1000))/obj.Attacker.FramePeriodicity_ms) + 1;
            attacker_frame_start_ms = obj.Attacker.FramePeriodicity_ms * (attacker_frame - 1) + obj.Attacker.timing_offset_us / 1000;

            %calculate the [t,f] values for the first attacker frame
            attacker_frame_end_ms = attacker_frame_start_ms + obj.Attacker.FramePeriodicity_ms * (attacker_frame);
            [attacker_tx,attacker_f_tx] = obj.Attacker.generate_frame_f_over_T_vals(attacker_frame_start_ms);
            [attacker_t_sampling,attacker_f_sampling] = obj.Attacker.generate_frame_f_over_Tsampling_vals(attacker_frame_start_ms);

            while attacker_frame_end_ms < victim_frame_end_ms
                % if the attacker frame ends before the victim frame, 
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

            %identify the points where the victim would actually observe
            %the intersection of the attacker and victim frames
            
           [attack_region_t, attack_region_f] = obj.identify_valid_frame_interference(victim_t_sampling, victim_f_sampling, attacker_tx, attacker_f_tx);

            %plot everything
            obj.plot_radar_frame(axis,obj.Victim,victim_tx,victim_f_tx,victim_t_sampling,victim_f_sampling);
            hold on
            obj.plot_radar_frame(axis,obj.Attacker,attacker_tx,attacker_f_tx,attacker_t_sampling,attacker_f_sampling);
            
            hold on
            p = scatter(axis,attack_region_t,attack_region_f,'.');
            p.DisplayName = sprintf('Valid attack region ');
            p.MarkerEdgeColor = 'green';
            hold off

            % by default, zoom the axis to only focus on the victim frame
            set(axis,'XLim',[victim_frame_start_ms,victim_frame_start_ms + (obj.Victim.ChirpCycleTime_us * 1e-3 * obj.Victim.NumChirps)]);
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
    
        function [attack_region_t,attack_region_f] = identify_valid_chirp_interference(obj, victim_t_sampling, victim_f_sampling, attacker_tx, attacker_f_tx)
            %{
                Purpose: identifies the points in the attacker chirp that
                    would actually interfere with the victim chirp (i.e: the
                    difference between the attacker and victim chirps at a
                    specific value of t is within half of the sampling
                    frequency of the victim
                Inputs:
                    [victim_t_sampling,victim_f_sampling]: the [t,f]
                        values of the sampling period for the victim
                    [attacker_tx, attacker_f_tx]: the [t,f] values of the
                        sampling period for the transmission period of the
                        attacker
                Outputs:
                    [attack_region_t, attack_region_f]: the [t,f] values
                    corresponding to the when the attacker is actually
                    interfering with the victim
            %}

            %identify intersection between victim_t_sampling and attacker_tx
            %appears to be a floating point discrepancy in the time values
            %and so try and reset them so that they are equal again

            victim_t_sampling = round(victim_t_sampling/obj.Victim.plotResolution_us) * obj.Victim.plotResolution_us;
            attacker_tx = round(attacker_tx / obj.Victim.plotResolution_us) * obj.Victim.plotResolution_us;

            [intersection_values,i_victim_t_sampling,i_attacker_tx] = intersect(victim_t_sampling,attacker_tx);
            intersection_victim_t = victim_t_sampling(i_victim_t_sampling);
            intersection_victim_f = victim_f_sampling(i_victim_t_sampling);
    
            intersection_attacker_t = attacker_tx(i_attacker_tx);
            intersection_attacker_f = attacker_f_tx(i_attacker_tx);
            
            %determine the maximum if frequency observable from the victim
            max_if_freq_GHz = (obj.Victim.ADC_SampleRate_MSps/2)/1000;
            valid_intersection = (abs(intersection_attacker_f - intersection_victim_f) < max_if_freq_GHz);
    
            attack_region_t = intersection_attacker_t(valid_intersection);
            attack_region_f = intersection_attacker_f(valid_intersection);
        end

        function [attack_region_t,attack_region_f] = identify_valid_frame_interference(obj, victim_t_sampling, victim_f_sampling, attacker_tx, attacker_f_tx)
            %{
                Purpose: identifies the points in the attacker chirp that
                    would actually interfere with the victim chirp (i.e: the
                    difference between the attacker and victim chirps at a
                    specific value of t is within half of the sampling
                    frequency of the victim
                Inputs:
                    [victim_t_sampling,victim_f_sampling]: the [t,f]
                        values of the sampling period for the victim
                    [attacker_tx, attacker_f_tx]: the [t,f] values of the
                        sampling period for the transmission period of the
                        attacker
                Outputs:
                    [attack_region_t, attack_region_f]: the [t,f] values
                    corresponding to the when the attacker is actually
                    interfering with the victim
            %}

            %identify intersection between victim_t_sampling and attacker_tx
            %appears to be a floating point discrepancy in the time values
            %and so try and reset them so that they are equal again
            
            %multiplying and dividing by 1000 to convert the resolution to
            %ms
            
            victim_t_sampling = round(victim_t_sampling*1000/obj.Victim.plotResolution_us) * obj.Victim.plotResolution_us/1000;
            attacker_tx = round(attacker_tx *1000 / obj.Victim.plotResolution_us) * obj.Victim.plotResolution_us/1000;

            [intersection_values,i_victim_t_sampling,i_attacker_tx] = intersect(victim_t_sampling,attacker_tx);
            intersection_victim_t = victim_t_sampling(i_victim_t_sampling);
            intersection_victim_f = victim_f_sampling(i_victim_t_sampling);
    
            intersection_attacker_t = attacker_tx(i_attacker_tx);
            intersection_attacker_f = attacker_f_tx(i_attacker_tx);
            
            %determine the maximum if frequency observable from the victim
            max_if_freq_GHz = (obj.Victim.ADC_SampleRate_MSps/2)/1000;
            valid_intersection = (abs(intersection_attacker_f - intersection_victim_f) < max_if_freq_GHz);
    
            attack_region_t = intersection_attacker_t(valid_intersection);
            attack_region_f = intersection_attacker_f(valid_intersection);
        end

    end
end