classdef FilterParameters < handle
    %FilterParameters: class containing the parameters of the state
    %estimation filter
    
    properties (Access=public)
        positionProcessNoiseStd
        attitudeProcessNoiseStd
        linVelProcessNoiseStd
        linAccProcessNoiseStd
        angVelProcessNoiseStd
        accBiasNoiseStd % Accelerometer bias noise standard deviation
        gyrBiasNoiseStd % Gyroscope bias noise standard deviation
        
        nIters % Number of iterations for the correction step (iterated EKF)
        temporalWindow % Temporal window of the state buffer [s]
        minCovMatrix
        maxCovMatrix
        minPosition
        maxPosition
        minAttitude
        maxAttitude
        minLinVel
        maxLinVel
        minAngVel
        maxAngVel
        minLinAcc
        maxLinAcc
        minAccBias
        maxAccBias
        minGyrBias
        maxGyrBias
        maxPositionStd
        maxAttitudeStd
        maxLinVelStd
        maxAngVelStd
        maxLinAccStd
        maxAccBiasStd
        maxGyrBiasStd
        gravityVector
    end
    
    methods (Access=public)
        function obj = FilterParameters(params)
            %FilterParameters Construct an instance of this class
            %   Detailed explanation goes here
            obj.positionProcessNoiseStd = params.positionProcessNoiseStd;
            obj.attitudeProcessNoiseStd = params.attitudeProcessNoiseStd;
            obj.linVelProcessNoiseStd = params.linVelProcessNoiseStd;
            obj.angVelProcessNoiseStd = params.angVelProcessNoiseStd;
            obj.linAccProcessNoiseStd = params.linAccProcessNoiseStd;
            obj.accBiasNoiseStd = params.accBiasNoiseStd;
            obj.gyrBiasNoiseStd = params.gyrBiasNoiseStd;
            % ...
        end

    end
end