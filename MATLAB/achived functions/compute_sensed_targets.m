
%only test with range or velocity at a time. If testing range, set velocity
%to be zero

function [estimated_ranges,estimated_velocities] = ...
    compute_sensed_targets(config_path,actual_range,actual_velocity,frames_to_compute,attack_start_frame)
    
    simulator = Simulator_revB();
    
    simulator.load_params_from_JSON(config_path);
    
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
    attack_type = "target";
    %attack_type = "target";
    %initialize the attacker
    simulator.Attacker.Subsystem_attacking.set_attack_mode(attack_type);
    
    %if it is desired to specify a specific attack location
    simulator.Attacker.Subsystem_attacking.set_desired_attack_location(actual_range,actual_velocity);

    simulator.run_simulation_attack_no_target(frames_to_compute,false);
    
    %get the return values

    estimated_ranges = simulator.Victim.Radar_Signal_Processor.range_estimates(attack_start_frame:frames_to_compute,1);
    estimated_velocities = simulator.Victim.Radar_Signal_Processor.velocity_estimates(attack_start_frame:frames_to_compute,1);
end

    