%%% ARCHIVED VERSION %%%%


classdef Target < handle
    

    properties (Access = public)
        distance_m
        speed_meters_per_s
        rcs_sq_meters
        operating_frequency_Hz

        %key FMCW parameters
        radar_target            %phased.RadarTarget object
        platform                %phased.Platform object

        %parameter to enable/disable the target
        enabled
    end

    properties (Access = private)
        propogation_speed

        %creating a position and velocity vector set for the target based
        %on the speed and distance parameters in the publicly available set
        position_m
        velocity_meters_per_s
    end

    methods
        function obj = Target(distance,speed_meters_per_s,rcs_sq_meters,operating_frequency_Hz, enabled)
            %{
                Purpose: creates an instance of the Target Class 
            %}
            obj.distance_m = distance;
            obj.rcs_sq_meters = rcs_sq_meters;
            obj.speed_meters_per_s = speed_meters_per_s;
            obj.operating_frequency_Hz = operating_frequency_Hz;
            obj.enabled = enabled;

            obj.propogation_speed = physconst('LightSpeed');
            
            
            %initialize the key FMCW parameters
            obj.configure_target_FMCW_params();

        end

        function  configure_target_FMCW_params(obj)
            %{
                Purpose: compute the relevant parameters needed for the
                    FMCW simulation
            %}
            obj.position_m = [obj.distance_m;0;0];
            obj.velocity_meters_per_s = [obj.speed_meters_per_s;0;0];

            obj.radar_target = phased.RadarTarget( ...
                'MeanRCS', obj.rcs_sq_meters, ...
                'PropagationSpeed',obj.propogation_speed, ...
                'OperatingFrequency', obj.operating_frequency_Hz);
            obj.platform = phased.Platform( ...
                'InitialPosition', obj.position_m, ...
                'Velocity',obj.velocity_meters_per_s);
        end
    end
end