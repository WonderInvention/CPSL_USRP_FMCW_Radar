%% For use in SIMULINK Model

%% Notes
%{
    -This class is originally based off of the Radar Class, but has been
    tailored to include things that are specific to the attacker
%}

classdef Attacker_revA < handle
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here

    properties (Access = public)
        %parameters for each of the subsystems in the class
        Subsystem_attacking         %a Subsystem_attacking class object
        Subsystem_tracking          %a Radar class object
        Subsystem_spectrum_sensing  %a Subsystem_spectrum_sensing class object
        
        %parameters for the position and speed of the attacker
        position_m
        velocity_m_per_s

        %parameter for the platform of the attacker
        platform            %phased.Platform object

    end

    methods (Access = public)
        function obj = Attacker_revA()
            %Radar construct an instance of this class
            %   Detailed explanation goes here
            obj.Subsystem_tracking = Radar_revA();
            obj.Subsystem_attacking = Subsystem_attacking();
            obj.Subsystem_spectrum_sensing = Subsystem_spectrum_sensing();
        end

        function configure_platform(obj)
            %{
                Prpose: configures the platform object as well as the
                platform objects for the subsystem modules as needed
            %}

            %for the attack module itself
            obj.platform = phased.Platform( ...
                'InitialPosition',obj.position_m, ...
                'Velocity',obj.velocity_m_per_s);

            %for the tracking subsystem
            obj.Subsystem_tracking.position_m = obj.position_m;
            obj.Subsystem_tracking.velocity_m_per_s = obj.velocity_m_per_s;
            obj.Subsystem_tracking.platform = obj.platform;
        end
    
        function set_desired_attack_parameters(obj, desired_range_m,desired_velocity_m_s)
            %{
                Purpose: sets the desired attack parameters for the
                attacker
            %}
            obj.Subsystem_attacking.desired_attack_range_m = desired_range_m;
            obj.Subsystem_attacking.desired_attack_velocity_m_s = desired_velocity_m_s;

        end
    
        function import_parameters_from_victim(obj, Victim)
        %{
            Purpose: This funciton allows the attacker to obtain its
            parameters from the victim. Basically, it allows me to use the
            attacker as a target emulator
        %}
            %imported parameters
        obj.Subsystem_attacking.StartFrequency_GHz = Victim.StartFrequency_GHz;
        obj.Subsystem_attacking.FrequencySlope_MHz_us= Victim.FrequencySlope_MHz_us;
        obj.Subsystem_attacking.IdleTime_us = Victim.IdleTime_us;                
        obj.Subsystem_attacking.RampEndTime_us = Victim.RampEndTime_us;                                  
        obj.Subsystem_attacking.Chirp_Tx_Bandwidth_MHz = Victim.Chirp_Tx_Bandwidth_MHz;          
        obj.Subsystem_attacking.ChirpCycleTime_us = Victim.ChirpCycleTime_us;
        obj.Subsystem_attacking.Lambda_m = Victim.Lambda_m;                        

        %frame parameters
        obj.Subsystem_attacking.NumChirps = Victim.NumChirps;
        obj.Subsystem_attacking.FramePeriodicity_ms = Victim.FramePeriodicity_ms;
        obj.Subsystem_attacking.ActiveFrameTime_ms = Victim.ActiveFrameTime_ms;              

        %parameters for the FMCW waveform, not setable parameters in UI
        obj.Subsystem_attacking.FMCW_sampling_rate_Hz = Victim.FMCW_sampling_rate_Hz;
        obj.Subsystem_attacking.sweep_time_s = Victim.sweep_time;
        obj.Subsystem_attacking.num_samples_per_frame = Victim.num_samples_per_frame;
        end
    
    end
end