classdef State < handle
    %State Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access=public)
        timestamp % Timestamp of the state (double) [s]
        stateVector % State vector (vector) [n x 1]
        covarianceMatrix % Associated covariance matrix (matrix) [n x n]
        measurements % Associated measurements (struct)
        corrected % If the state has been corrected with the latest measurement (bool)
    end
    
    methods (Access=public)
        function obj = State(ts, x, P)
            %State Construct an instance of this class
            %   Detailed explanation goes here
            
            assert(ts >= 0.0, 'State::State: timestamp must be >= 0!');
            assert(length(x)-1 == size(P, 1), 'State::State: x and P dimension mismatch!');
            assert(size(P, 1) == size(P, 2), 'State::State: P must be square!');
            
            obj.timestamp = ts;
            obj.stateVector = x;
            obj.covarianceMatrix = P;
            obj.measurements = {};
            obj.corrected = true;
        end
        
        function addMeasurement(obj, z)
            %addMeasurement: add a measurement to the state
            nMeasurements = length(obj.measurements); % Current nb of measurements
            
            obj.measurements{nMeasurements + 1} = z;
            obj.corrected = false; % Reset corrected flag!
        end
        
        function n = getNbMeasurements(obj)
            % getNbMeasurements: return the number of measurements associated to the state
            n = length(obj.measurements);
        end
        
        function enforceCovMatSymmetry(obj)
            % enforceCovMatSymmetry
            obj.covarianceMatrix = 0.5 * (obj.covarianceMatrix + obj.covarianceMatrix');
        end
        
        function isGood = checkStateVector(obj)
            % checkStateVector: check the current state vector
            
            isGood = true;
            % Check for NaN
            maskNan = isnan(obj.stateVector);
            idxs = find(maskNan == 1);
            if ~isempty(idxs)
                %warning('State::checkStateVector: isnan detected!');
                obj.stateVector(maskNan) = 0.0;
                isGood = false;
            end
            
            % Check for inf
            maskInf = isinf(obj.stateVector);
            idxs = find(maskInf == 1);
            if ~isempty(idxs)
                %warning('State::checkStateVector: isinf detected!');
                obj.stateVector(maskInf) = 0.0;
                isGood = false;
            end
        end
        
        function isGood = checkCovarianceMatrix(obj)
            % checkCovarianceMatrix: check the current covariance matrix
            
            isGood = true;
            % Check for NaN
            maskNan = isnan(obj.covarianceMatrix);
            idxs = find(maskNan == 1);
            if ~isempty(idxs)
                %warning('State::covarianceMatrix: isnan detected!');
                obj.covarianceMatrix(maskNan) = 0.0;
                isGood = false;
            end
            
            % Check for inf
            maskInf = isinf(obj.covarianceMatrix);
            idxs = find(maskInf == 1);
            if ~isempty(idxs)
                %warning('State::covarianceMatrix: isinf detected!');
                obj.covarianceMatrix(maskInf) = 0.0;
                isGood = false;
            end
        end
        
        function str = toString(obj)
            % toString: convert object to string for debug information
            str = ['State @ t=', num2str(obj.timestamp), ' s with ', ...
                num2str(length(obj.measurements)), ' measurements, x=(', ...
                num2str(obj.stateVector'), ')'];
        end
        
        % Overload comparison operators
        function res = lt(obj1, obj2)
            if obj1.timestamp < obj2.timestamp
                res = true;
            else
                res = false;
            end
        end
        
        function res = gt(obj1, obj2)
            if obj1.timestamp > obj2.timestamp
                res = true;
            else
                res = false;
            end
        end
        
        function res = le(obj1, obj2)
            if obj1.timestamp <= obj2.timestamp
                res = true;
            else
                res = false;
            end
        end
        
        function res = ge(obj1, obj2)
            if obj1.timestamp >= obj2.timestamp
                res = true;
            else
                res = false;
            end
        end
        
        function res = eq(obj1, obj2)
            if (obj1.timestamp == obj2.timestamp) ...
                    && strcmp(obj1.sensorName, obj2.sensorName)
                res = true;
            else
                res = false;
            end
        end
        
        function res = ne(obj1, obj2)
            res = ~(obj1 == obj2);
        end
    end
end

