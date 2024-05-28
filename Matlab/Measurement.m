classdef Measurement < handle
    %Measurement Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access=public)
        timestamp % Timestamp of the measurement [s] (double)
        measurementVector % Measurement vector (vector) [nx1]
        covarianceMatrix % Measurement covariance matrix (matrix) [nxn]
        id % Measurement unique ID (uint64)
        sensorName % Name of the corresponding sensor (string)
        mahalanobisThresh % Mahalanobis threshold (double)
        %sensorTransformation % Transformation of the sensor wrt to origin (matrix [4x4])
    end
    
    methods (Access=public)
        function obj = Measurement(ts, z, R, id, sensorName, mahaThresh)
            %State Construct an instance of this class
            %   Detailed explanation goes here
            assert(ts >= 0.0, 'Measurement::Measurement: timestamp must be >= 0!');
            assert(length(z) == size(R, 1), 'Measurement::Measurement: z and R dimension mismatch!');
            assert(size(R, 1) == size(R, 2), 'Measurement::Measurement: R must be square!');
            assert(mahaThresh > 0.0, 'Measurement::Measurement: mahaThresh must be > 0!');
            
            obj.timestamp = ts;
            obj.measurementVector = z;
            obj.covarianceMatrix = R;
            obj.id = id;
            obj.sensorName = sensorName;
            obj.mahalanobisThresh = mahaThresh;
        end
        
        function str = toString(obj)
            %toString Summary of this method goes here
            %   Detailed explanation goes here
            str = ['Measurement from ',  obj.sensorName, ': id=', ...
                num2str(obj.id), ', @ t=', num2str(obj.timestamp), ...
                ' s, z=(', num2str(obj.measurementVector'), ...
                '), R=(', num2str(diag(obj.covarianceMatrix)'), ')'];
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