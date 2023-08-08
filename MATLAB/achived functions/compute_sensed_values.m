

function [actual_slope,estimated_slope,...
    actual_chirp_duration,estimated_chirp_duration,...
    actual_frame_duration, estimated_frame_duration] = ...
    compute_sensed_values(config_path,slope_MHz_us,chirp_cycle_period_us, frames_to_compute)
    
    simulator = Simulator_revB();
    
    simulator.load_params_from_JSON(config_path);
    
    simulator.Victim.FrequencySlope_MHz_us = slope_MHz_us;
    simulator.Victim.ChirpCycleTime_us = chirp_cycle_period_us;
    simulator.Victim.compute_calculated_vals();

    %set attacker parameters
    simulator.Attacker.Subsystem_tracking.FrequencySlope_MHz_us = slope_MHz_us;
    simulator.Attacker.Subsystem_tracking.ChirpCycleTime_us = chirp_cycle_period_us;
    simulator.Attacker.Subsystem_tracking.compute_calculated_vals();
    
    %apply timing offsets as desired
    simulator.Victim.timing_offset_us = 0;
    simulator.Attacker.Subsystem_tracking.timing_offset_us = 0;
    
    %configure the FMCW parameters
    simulator.configure_FMCW_Radar_parameters();
    
    %load default attacker, and victim positions and velocities
    simulator.load_realistic_attacker_and_victim_position_and_velocity();
    
    %load the target
    simulator.load_target_realistic(50,15);

    %initialize victim and simulation parameters

    
    %specify whether or not to record a movie of the range-doppler plot
    record_movie = false;
    simulator.Victim.Radar_Signal_Processor.configure_movie_capture(frames_to_compute,record_movie);
    
    %pre-compute the victim's chirps
    simulator.Victim.precompute_radar_chirps();

    %initialize the attacker parameters
    simulator.Attacker.initialize_attacker(...
                            simulator.Victim.FMCW_sampling_rate_Hz * 1e-6,...
                            simulator.Victim.StartFrequency_GHz,...
                            simulator.Victim.Chirp_Tx_Bandwidth_MHz);
    
    %initialize the sensing subsystem's debugger
    simulator.Attacker.Subsystem_spectrum_sensing.initialize_debugger(0,simulator.Victim,frames_to_compute);
    %set to change 1 to zero to disable
    
    %specify the type of emulation ("target", 
            % "velocity spoof - noisy", 
            % "velocity spoof - similar velocity",
            % "range spoof - similar slope")
    
    %attack_type = "range spoof - similar slope";
    attack_type = "range spoof - similar slope,velocity spoof - noisy";
    %attack_type = "target";
    %initialize the attacker
    simulator.Attacker.Subsystem_attacking.set_attack_mode(attack_type);
    
    %if it is desired to specify a specific attack location
    %simulator.Attacker.Subsystem_attacking.set_desired_attack_location(100,7);

    simulator.run_simulation_with_attack(frames_to_compute,false);
    
    %get the return values

    %frame duration
    estimated_frame_duration = simulator.Attacker.Subsystem_spectrum_sensing.frame_tracking.average_frame_duration * 1e-3;
    actual_frame_duration = simulator.Victim.FramePeriodicity_ms;
    %chirp duration
    estimated_chirp_duration = simulator.Attacker.Subsystem_spectrum_sensing.frame_tracking.average_chirp_duration;
    actual_chirp_duration = simulator.Victim.ChirpCycleTime_us;
    %slope
    estimated_slope = simulator.Attacker.Subsystem_spectrum_sensing.frame_tracking.average_slope;
    actual_slope = simulator.Victim.FrequencySlope_MHz_us;
end

    