classdef Ekf < handle
    %Ekf: Extended Kalman Filter for state estimation
    % the state vector (of length 22) is composed with:
    % - 3D position (in ENU local frame) [m]
    % - 3D attitude (quaternion) (from body to ENU local frame)
    % - 3D linear velocity (in body frame) [m/s]
    % - 3D angular velocity (in body frame) [rad/s]
    % - 3D linear acceleration (in body frame) [m/s2]
    % - 3D accelerometer bias (in sensor frame) [m/s2]
    % - 3D gyroscope bias (in sensor frame) [rad/s]
    % - Optional: 3D magnetometer bias (in sensor frame) [uT]
    
    properties (Access=public)
        %states % Array/struct of states
        stateVector % State vector (size [22x1])
        stateVectorPred % State vector after prediction (and before correction) (size [22x1])
        covarianceMatrix % Covariance matrix (size [21x21])
        covarianceMatrixPred % Covariance matrix after prediction (and before correction)  (size [21x21])
        F % State transition matrix (size [21x21])
        W % Process noise Jacobian matrix
        Q % Process noise matrix
        G % Matrix for resetting the covariance matrix (size [21x21])
        I % Identity matrix (used during correction steps) (size [21x21])
        isInit % Boolean flag indicating if the filter is initialized
        isSane % Boolean flag indicating if the filter estimation is valid
        stateSanityMask % Bitfield mask to track which state is a problem in the state vector
        covSanityMask % Bitfield mask to track which state is a problem in the covariance matrix
        % Parameters
        %parameters % Configuration parameters of the EKF
        % Noises
        positionProcessNoiseStd % Process/prediction noise on the position
        attitudeProcessNoiseStd % Process/prediction noise on the attitude
        linVelProcessNoiseStd % Process/prediction noise on the linear velocity
        angVelProcessNoiseStd % Process/prediction noise on the angular velocity
        linAccProcessNoiseStd % Process/prediction noise on the linear acceleration
        accBiasNoiseStd % Process/prediction noise on the accelero bias
        gyrBiasNoiseStd % Process/prediction noise on the gyro bias

        % Covariance matrix bounds
        minCovariance
        maxCovariance
        maxPositionStd
        maxAttitudeStd
        maxLinVelStd
        maxAngVelStd
        maxLinAccStd
        maxAccBiasStd
        maxGyrBiasStd
        
        % State vector bounds
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
        
        % Vehicle parameters
        wheelRadius % Vehicle wheel radius [m]
        distRearAxleToCog % Distance between rear axle and vehicle CoG [m]
        distFrontAxleToCog % Distance between front axle and vehicle CoG [m]
        
        % Others
        gravity % Earth gravity vector [m/s2]
        nStates
        nErrorStates
    end
    
    methods (Access=public)
        function obj = Ekf()
            %Ekf Construct an instance of this class
            
            deg2rad = pi / 180.0;

            % Set parameters
            obj.nStates = 22; % Number of states
            obj.nErrorStates = 21; % Number of error states
            
            % Set process/prediction noise
            obj.positionProcessNoiseStd = [1e-4 / 3.; 1e-4 / 3.; 1e-4 / 3.]; % [m]
            obj.attitudeProcessNoiseStd = [0.1*deg2rad/3.0; 0.1*deg2rad/3.0; 0.1*deg2rad/3.0]; % [rad/s]
            obj.linVelProcessNoiseStd = [0.1/3.0; 0.1/3.0; 0.1/3.0]; % [m/s]
            obj.angVelProcessNoiseStd = [0.1*deg2rad/3.0; 0.1*deg2rad/3.0; 0.1*deg2rad/3.0]; % [rad/s]
            obj.linAccProcessNoiseStd = [1e-2/3.0; 1e-2/3.0; 1e-2/3.0]; % [m/s2]
            obj.accBiasNoiseStd = [1e-4/3.0; 1e-4/3.0; 1e-4/3.0]; % [m/s2]
            obj.gyrBiasNoiseStd = [1e-5/3.0; 1e-5/3.0; 1e-5/3.0]; % [rad/s]
            
            % Set state vector bounds
            obj.minPosition = [-10e3, -10e3, -10e3]; % [m]
            obj.maxPosition = [10e3, 10e3, 10e3]; % [m]
            obj.minAttitude = [-pi; -pi; -pi]; % [rad]
            obj.maxAttitude = [pi; pi; pi]; % [rad]
            obj.minLinVel = [-30., -10.0, -10.0]; % [m/s]
            obj.maxLinVel = [30., 10.0, 10.0]; % [m/s]
            obj.minAngVel = [-7.0, -7.0, -7.0]; % [rad/s]
            obj.maxAngVel = [7.0, 7.0, 7.0]; % [rad/s]
            obj.minLinAcc = [-10.0, -4.0, -4.0]; % [m/s2]
            obj.maxLinAcc = [8.0, 4.0, 4.0]; % [m/s2]
            obj.minAccBias = [-10.0, -4.0, -4.0]; % [m/s2]
            obj.maxAccBias = [8.0, 4.0, 4.0]; % [m/s2]
            obj.minGyrBias = [-7.0, -7.0, -7.0]; % [rad/s]
            obj.maxGyrBias = [7.0, 7.0, 7.0]; % [rad/s]
            
            % Set covariance matrix bounds
            obj.minCovariance = 1e-5;
            obj.maxCovariance = 1e3;
            
            obj.maxPositionStd = [10.0, 10.0, 10.0]; % [m]
            obj.maxAttitudeStd = [0.5, 0.5, 0.5]; % [rad]
            obj.maxLinVelStd = [1.0, 1.0, 1.0]; % [m/s]
            obj.maxAngVelStd = [0.2, 0.2, 0.2]; % [rad/s]
            obj.maxLinAccStd = [0.5, 0.5, 0.5]; % [m/s2]
            obj.maxAccBiasStd = [0.3, 0.3, 0.3]; % [m/s2]
            obj.maxGyrBiasStd = [0.1, 0.1, 0.1]; % [rad/s]
            
            % Set parameters
            obj.wheelRadius = 0.3575; % [m]
            obj.distRearAxleToCog = 1.135; % [m]
            obj.distFrontAxleToCog = 1.115; % [m]
            
            obj.gravity = [0.0; 0.0; -9.81]; % Earth gravity vector [m/s2]
            
            obj.stateVector = zeros(22, 1); % Init state
            obj.stateVector(4) = 1.0; % For quaternion qw !!!
            obj.stateVectorPred = zeros(22, 1); % Init state
            obj.stateVectorPred(4) = 1.0; % For quaternion qw !!!
            obj.covarianceMatrix = eye(21, 21);
            obj.covarianceMatrixPred = eye(21, 21);
            obj.F = eye(21, 21);
            obj.Q = zeros(12, 12);
            obj.W = zeros(21, 12);
            obj.G = eye(21, 21);
            obj.I = eye(21, 21); % Identity matrix
            
            obj.isInit = false;
            obj.isSane = false;
            obj.stateSanityMask = zeros(7, 1);
            obj.covSanityMask = zeros(7, 1);
        end
        
        function init(obj, ts0, x0, P0)
            % init: Initialize the Kalman filter with first state and covariance
            % Args:
            % - ts0 (double): timestamp of the state [s]
            % - x0 (vector): initial state vector
            % - P0 (matrix): initial covariance matrix
            
            assert(ts0 >= 0.0, 'Ekf::init: timestamp must be >= 0!');
            assert(length(x0) == obj.nStates, 'Ekf::init: wrong size for x0!')
            assert(size(P0, 1) == obj.nErrorStates, 'Ekf::init: wrong size for P0!');
            assert(size(P0, 1) == size(P0, 2), 'Ekf::init: P0 must be a square matrix!');
            
            assert(~any(isnan(x0), 'all'), 'Ekf::init: NaN in x0');
            assert(~any(isinf(x0), 'all'), 'Ekf::init: INF in x0');
            
            assert(~any(isnan(P0), 'all'), 'Ekf::init: NaN in P0');
            assert(~any(isinf(P0), 'all'), 'Ekf::init: INF in P0');
            
            obj.stateVector = x0; % Reset states
            obj.covarianceMatrix = P0;
            
            obj.isInit = true;
            obj.isSane = true; % Assume sane
            obj.stateSanityMask = ones(7, 1);
            obj.covSanityMask = ones(7, 1);
        end
        
        function reset(obj)
            % reset: reset the EKF filter
            
            %warning('Ekf::reset: resetting the EKF!');
            %fprintf('Ekf::reset: resetting the EKF!\n');
            
            obj.stateVector = zeros(22, 1); % Init state
            obj.stateVector(4) = 1.0; % For quaternion qw !!!
            obj.stateVectorPred = zeros(22, 1); % Init state
            obj.stateVectorPred(4) = 1.0; % For quaternion qw !!!
            obj.covarianceMatrix = eye(21, 21);
            obj.covarianceMatrixPred = eye(21, 21);
            
            obj.isInit = false;
            obj.isSane = false;
            
        end
        
        function success = predict(obj, dt)
            % predict: EKF prediction step, predict next state given the
            % current state and system dynamics
            % Args:
            % - dt (double): delta time, period of integration [s]
            % Return:
            % - success (bool): boolean flag indicating the prediction was successful
            
            success = false;
            
            if ~obj.isInit % Check if filter is initialized
                %warning('Ekf::predict: filter not initialized!');
                %fprintf('Ekf::predict: filter not initialized!\n');
                obj.stateVectorPred = obj.stateVector;
                obj.covarianceMatrixPred = obj.covarianceMatrix;
                return;
            end
            
            % Sampling interval [s]
            assert(dt > 0.0, 'Ekf::predict: dt <= 0!');
            
            % Get current state
            pCurr = obj.stateVector(1:3); % Get current position
            qCurr = obj.stateVector(4:7); % Get current quaternion
            vCurr = obj.stateVector(8:10); % Get current lin vel
            wCurr = obj.stateVector(11:13); % Get current ang vel
            aCurr = obj.stateVector(14:16); % Get current lin acc
            baCurr = obj.stateVector(17:19); % Get current acc bias
            bgCurr = obj.stateVector(20:22); % Get current gyr bias
            
            % --- Predict next state ---
            
            % Position: p_{k+1} = p_{k} + R_k v_k dt + 0.5 R_k a_k dt^2
            R = quatToRotMat(qCurr); % Transform quaternion to rotation matrix
            obj.stateVectorPred(1:3) = pCurr + R * (vCurr * dt + 0.5 * aCurr * dt^2);
            
            % Attitude: q_{k+1} = q_{k} @ q(w * dt)
            deltaAngle = wCurr * dt; % w *dt [rad]
            qRot = deltaAngleToQuaternion(deltaAngle);
            qNext = quaternionProduct(qCurr, qRot);
            qNext = quaternionNormalize(qNext); % Normalize quaternion to keep a unit quaternion
            obj.stateVectorPred(4:7) = qNext(1:4);
            
            % Linear velocity: v_{k+1} = v_{k} + a_{k} dt
            obj.stateVectorPred(8:10) = vCurr + aCurr * dt;
            
            % Angular velocity: w_{k+1} = w_{k}
            obj.stateVectorPred(11:13) = wCurr;
            
            % Linear acceleration: a_{k+1} = a_{k}
            obj.stateVectorPred(14:16) = aCurr;
            
            % Accelero bias: ba_{k+1} = ba_{k}
            obj.stateVectorPred(17:19) = baCurr;
            
            % Gyro bias: bg_{k+1} = bg_{k}
            obj.stateVectorPred(20:22) = bgCurr;
            
            % --- State transition matrix ---
            % 1:3 -> position, 4:6 -> attitude, 7:9 -> lin vel, 
            % 10:12 -> ang vel, 13:15 -> lin acc, 16:18 -> acc bias,
            % 19:21 -> gyr bias
            Vx = skewMatrix(wCurr);
            
            %obj.F = eye(21, 21);
            % dp/dp
            %obj.F(1:3, 1:3) = eye(3);
            % dp/dt
            obj.F(1:3, 4:6) = -R * Vx * dt;
            % dp/dv
            obj.F(1:3, 7:9) = R * dt;
            
            % dt/dt
            obj.F(4:6, 4:6) = quatToRotMat(qRot);
            % dt/dw
            obj.F(4:6, 10:12) = dt * eye(3);
            
            % dv/dv
            %obj.F(7:9, 7:9) = eye(3);
            % dv/da
            obj.F(7:9, 13:15) = dt * eye(3);
            
            % dw/dw
            %obj.F(10:12, 10:12) = eye(3);
            
            % da/da
            %obj.F(13:15, 13:15) = eye(3);
            
            % dba/dba
            %obj.F(16:18, 16:18) = eye(3);
            
            % dbg/dbg
            %obj.F(19:21, 19:21) = eye(3);
            
            assert(~any(isnan(obj.F), 'all'), 'Ekf::predict: NaN detected in F!');
            assert(~any(isinf(obj.F), 'all'), 'Ekf::predict: INF detected in F!');
            
            % --- Process noise covariance matrix ---
            % --- Propagate covariance matrix: P = F * P * F^t + W * Q * W^t ---
            
            % V1
            %obj.Q = zeros(21, 21);
            %obj.Q(1:3, 1:3) = eye(3) * obj.positionProcessNoiseStd^2;
            %obj.Q(4:6, 4:6) = eye(3) * obj.attitudeProcessNoiseStd^2;
            %obj.Q(7:9, 7:9) = eye(3) * obj.linVelProcessNoiseStd^2;
            %...
            %obj.covarianceMatrixPred = obj.F * obj.covarianceMatrix * obj.F' + obj.Q;
            
            % V2
            obj.covarianceMatrixPred = obj.F * obj.covarianceMatrix * obj.F';
            obj.covarianceMatrixPred(7:9, 7:9) = ...
                obj.covarianceMatrixPred(7:9, 7:9) + diag(obj.linVelProcessNoiseStd.^2);
            obj.covarianceMatrixPred(10:12, 10:12) = ...
                obj.covarianceMatrixPred(10:12, 10:12) + diag(obj.angVelProcessNoiseStd.^2);
            obj.covarianceMatrixPred(13:15, 13:15) = ...
                obj.covarianceMatrixPred(13:15, 13:15) + diag(obj.linAccProcessNoiseStd.^2);
            obj.covarianceMatrixPred(16:18, 16:18) = ...
                obj.covarianceMatrixPred(16:18, 16:18) + diag(obj.accBiasNoiseStd.^2);
            obj.covarianceMatrixPred(19:21, 19:21) = ...
                obj.covarianceMatrixPred(19:21, 19:21) + diag(obj.gyrBiasNoiseStd.^2);
            
            obj.stateVector = obj.stateVectorPred; % Copy vector
            obj.covarianceMatrix = obj.covarianceMatrixPred; % Copy matrix
            
            % Check state vector and covariance matrix
            success = obj.checkStateVector() && obj.checkCovarianceMatrix();
        end
        
        function success = correct(obj, z, R, mahaThresh, z_pred, H, sensor)
            % correct: generic correction step function
            % Args:
            % - z (vector): measurement vector
            % - R (matrix): measurement noise covariance matrix
            % - mahaThresh (double): Mahalanobis distance threshold
            % - z_pred (vector): predicted measurement vector
            % - H (matrix): measurement Jacobian matrix
            % - sensor (string): name of the sensor
            % Return:
            % - success (bool): boolean flag indicating if the correction was successful
            
            assert(length(z) == length(z_pred), 'Ekf::correct (%s): len(z) != len(z_pred)', sensor);
            assert(length(z) == size(R, 1) && length(z) == size(R, 2), 'Ekf::correct (%s): size(R) error!', sensor);
            assert(length(z) == size(H, 1) && size(obj.covarianceMatrix, 1) == size(H, 2), 'Ekf::correct (%s): size(H) error!', sensor);
            
            success = false;
            
            % Get covariance matrix
            P = obj.covarianceMatrix;

            % Innovation vector
            nu = z - z_pred;
            
            % Innovation matrix
            PHt = P * H';
            S = H * PHt + R;
            S_inv = inv(S);
            
            % Mahalanobis outlier test
            maha = nu' * S_inv * nu;
            if (maha < 0.0)
                str = sprintf('Ekf::correct: maha < 0.0 (%s)!\n', sensor);
                %disp(str);
                fprintf(str);
                return;
            elseif (maha > mahaThresh)
                % Check if we can still continue applying the measurement but with
                % much higher covariance (This is to avoid divergence of the solution).
                if (maha <= 3.0 * mahaThresh)
                    str = sprintf('Ekf::correct: Mahalanobis test failed: attenuating confidence (%s)!\n', sensor);
                    %disp(str);
                    fprintf(str);
                    S = H * PHt + R * (maha / mahaThresh);
                    S_inv = inv(S);
                else
                    str = sprintf('Ekf::correct: Mahalanobis test failed (%s)!\n', sensor);
                    %disp(str);
                    fprintf(str);
                    return;
                end
            end
            
            % Kalman gain
            K = PHt * S_inv;
            
            % State vector correction
            dx = K * nu; % Error-state vector
            assert(~any(isnan(dx)), 'Ekf::correct (%s): NAN detected in dx!', sensor);
            assert(~any(isinf(dx)), 'Ekf::correct (%s): INF detected in dx!', sensor);
            obj.stateVector = obj.addErrorState(obj.stateVector, dx);
                        
            % Covariance matrix correction
            IKH = obj.I - K * H;
            P = IKH * P * IKH' + K * R * K';
            obj.covarianceMatrix = obj.errorReset(P, dx);
            
            % Check state vector and covariance matrix
            success = obj.checkStateVector() && obj.checkCovarianceMatrix();
        end
        
        function success = correctImu(obj, z, R, mahaThresh)
            % correctImu: EKF correction step with IMU measurement
            % Args:
            % - z (vector): measurement vector (first 3x1 accelero then 3x1 gyro)
            % - R (matrix): measurement noise covariance matrix
            % - mahaThresh (double): Mahalanobis distance threshold
            % Return:
            % - success (bool): boolean flag indicating if the correction was successful
            
            assert(length(z) == 6, 'Ekf::correctImu: z must be of size [6x1]!');
            assert(size(R, 1) == 6 && size(R, 2) == 6, 'Ekf::correctImu: R must be of size [6x6]!');
            assert(mahaThresh > 0.0, 'Ekf::correctImu: mahaThresh must be > 0!');

            % Predict measurement from current state
            q = obj.stateVector(4:7);
            Rot = quatToRotMat(q);
            Rtg = Rot' * obj.gravity;
            z_pred = zeros(6, 1);
            z_pred(1:3, 1) = obj.stateVector(14:16) + Rtg + obj.stateVector(17:19);
            z_pred(4:6, 1) = obj.stateVector(11:13) + obj.stateVector(20:22);

            % Measurement function Jacobian matrix
            H = zeros(6, 21);
            H(1:3, 4:6) = skewMatrix(Rtg);
            H(1:3, 13:15) = eye(3);
            H(1:3, 16:18) = eye(3);
            H(4:6, 10:12) = eye(3);
            H(4:6, 19:21) = eye(3);
            
            success = obj.correct(z, R, mahaThresh, z_pred, H, 'IMU');
        end
        
        function correctPosition(obj, z, R, mahaThresh)
            % correctPosition: EKF correction step with position measurement (eg GPS)
            % Args:
            % - z (vector): measurement vector
            % - R (matrix): measurement noise covariance matrix
            % - mahaThresh (double): Mahalanobis distance threshold
            % Return:
            % - success (bool): boolean flag indicating if the correction was successful
            
            assert(length(z) == 3, 'Ekf::correctPosition: z must be of size [6x1]!');
            assert(size(R, 1) == 3 && size(R, 2) == 3, 'Ekf::correctPosition: R must be of size [6x6]!');
            assert(mahaThresh > 0.0, 'Ekf::correctPosition: mahaThresh must be > 0!');
            
            % Predict measurement from current state
            z_pred = obj.stateVector(1:3);

            % Measurement function Jacobian matrix
            H = zeros(3, 21);
            H(1:3, 1:3) = eye(3);
            
            success = obj.correct(z, R, mahaThresh, z_pred, H, 'GPS POS');
        end
        
        function correctLinAngVel(obj, z, R, mahaThresh)
            % correctLinAngVel: EKF correction step with linear and angular velocity measurements
            % Args:
            % - z (vector): measurement vector
            % - R (matrix): measurement noise covariance matrix
            % - mahaThresh (double): Mahalanobis distance threshold
            % Return:
            % - success (bool): boolean flag indicating if the correction was successful
            
            assert(length(z) == 6, 'Ekf::correctLinAngVel: z must be of size [6x1]!');
            assert(size(R, 1) == 6 && size(R, 2) == 6, 'Ekf::correctLinAngVel: R must be of size [6x6]!');
            assert(mahaThresh > 0.0, 'Ekf::correctLinAngVel: mahaThresh must be > 0!');
            
            % Predict measurement from current state
            z_pred = zeros(6, 1);
            z_pred(1:3) = obj.stateVector(8:10); % Linear velocity
            z_pred(4:6) = obj.stateVector(11:13); % Angular velocity

            % Measurement function Jacobian matrix
            H = zeros(6, 21);
            H(1:3, 7:9) = eye(3);
            H(4:6, 10:12) = eye(3);
            
            success = obj.correct(z, R, mahaThresh, z_pred, H, 'LinAngVel');
        end
        
        function correctLinAcc(obj, z, R, mahaThresh)
            % correctLinAcc: EKF correction step with linear acceleration measurement
            % Args:
            % - z (vector): measurement vector
            % - R (matrix): measurement noise covariance matrix
            % - mahaThresh (double): Mahalanobis distance threshold
            % Return:
            % - success (bool): boolean flag indicating if the correction was successful
            
            assert(length(z) == 3, 'Ekf::correctLinAcc: z must be of size [6x1]!');
            assert(size(R, 1) == 3 && size(R, 2) == 3, 'Ekf::correctLinAcc: R must be of size [6x6]!');
            assert(mahaThresh > 0.0, 'Ekf::correctLinAcc: mahaThresh must be > 0!');
            
            % Predict measurement from current state
            z_pred = obj.stateVector(14:16);

            % Measurement function Jacobian matrix
            H = zeros(3, 21);
            H(1:3, 13:15) = eye(3);
            
            success = obj.correct(z, R, mahaThresh, z_pred, H, 'LinAcc');
        end
        
        function correctYaw(obj, z, R, mahaThresh)
            % correctYaw: EKF correction step with yaw/heading angle measurement
            % Args:
            % - z (vector): measurement vector
            % - R (matrix): measurement noise covariance matrix
            % - mahaThresh (double): Mahalanobis distance threshold
            % Return:
            % - success (bool): boolean flag indicating if the correction was successful
            
            assert(length(z) == 1, 'Ekf::correctYaw: z must be of size [1x1]!');
            assert(size(R, 1) == 1 && size(R, 2) == 1, 'Ekf::correctYaw: R must be of size [1x1]!');
            assert(mahaThresh > 0.0, 'Ekf::correctYaw: mahaThresh must be > 0!');
            
            assert(abs(z(1)) <= pi, 'Ekf::correctYaw: z=yaw must be in [-pi, pi]!');
            
            % Predict measurement from current state
            q = obj.stateVector(4:7);
            %rpy = quatToRollPitchYaw(q);
            %z_pred = rpy(3);
            z_pred = atan2(2.0*(q(1)*q(4) + q(2)*q(3)), 1.0 - 2.0*(q(3)^2 + q(4)^2));

            % Measurement function Jacobian matrix
            % Hx [1x4]
            Hx = [-(2*q(4))/(((4*(q(1)*q(4) + q(2)*q(3))^2)/(2*q(3)^2 + 2*q(4)^2 - 1)^2 + 1)*(2*q(3)^2 + 2*q(4)^2 - 1)), ...
                -(2*q(3))/(((4*(q(1)*q(4) + q(2)*q(3))^2)/(2*q(3)^2 + 2*q(4)^2 - 1)^2 + 1)*(2*q(3)^2 + 2*q(4)^2 - 1)), ...
                -((2*q(2))/(2*q(3)^2 + 2*q(4)^2 - 1) - (4*q(3)*(2*q(1)*q(4) + 2*q(2)*q(3)))/(2*q(3)^2 + 2*q(4)^2 - 1)^2)/((2*q(1)*q(4) + 2*q(2)*q(3))^2/(2*q(3)^2 + 2*q(4)^2 - 1)^2 + 1), ...
                -((2*q(1))/(2*q(3)^2 + 2*q(4)^2 - 1) - (4*q(4)*(2*q(1)*q(4) + 2*q(2)*q(3)))/(2*q(3)^2 + 2*q(4)^2 - 1)^2)/((2*q(1)*q(4) + 2*q(2)*q(3))^2/(2*q(3)^2 + 2*q(4)^2 - 1)^2 + 1)];
            % Qdt [4x3]
            Qdt = 0.5 * [-q(2), -q(3), -q(4);
                q(1), -q(4), q(3);
                q(4), q(1), -q(2);
                -q(3), q(2), q(1)];
            H = zeros(1, 21);
            H(1, 4:6) = Hx * Qdt;
            
            success = obj.correct(z, R, mahaThresh, z_pred, H, 'YAW');
        end
        
        function x = addErrorState(obj, x, dx)
            % addErrorState: add the error state correction to the state vector
            % Args:
            % - x (vector): state vector to correct
            % - dx (vector): error-state vector
            % Return:
            % - x (vector): corrected state vector
            
            % 3D position
            x(1:3) = x(1:3) + dx(1:3);
            
            % 3D attitude
            qCurr = x(4:7);
            deltaAngle = dx(4:6);
            qRot = deltaAngleToQuaternion(deltaAngle);
            qNew = quaternionProduct(qCurr, qRot);
            qNew = quaternionNormalize(qNew);
            x(4:7) = qNew(1:4);
            
            % 3D linear velocity
            x(8:10) = x(8:10) + dx(7:9);
            
            % 3D angular velocity
            x(11:13) = x(11:13) + dx(10:12);
            
            % 3D linear acceleration
            x(14:16) = x(14:16) + dx(13:15);
            
            % Accelero bias
            x(17:19) = x(17:19) + dx(16:18);
            
            % Gyro bias
            x(20:22) = x(20:22) + dx(19:21);
        end
        
        function P = errorReset(obj, P, dx)
            % errorReset: reset the error in the covariance matrix
            % Args:
            % - P (matrix): covariance matrix to reset
            % - dx (vector): error-state vector
            % Return:
            % - P (matrix): resetted covariance matrix
            
            attError = dx(4:6);
            if norm(attError) < 1e-9 % No attitude error, don't need to reset covariance matrix
                return;
            end
            
            obj.G(4:6, 4:6) = eye(3) - 0.5 * skewMatrix(attError);
            
            P = obj.G * P * obj.G';
        end
        
        function isGood = checkSanity(obj)
            % checkSanity: check the filter sanity (state vector and covariance matrix)
            
            % --- Check state vector ---
            
            % Reset mask
            obj.stateSanityMask = zeros(7, 1);
            
            x = obj.stateVector;
            
            % Check 3D position
            if x(1) > obj.maxPosition(1) || x(1) < obj.minPosition(1) ...
                    || x(2) > obj.maxPosition(2) || x(2) < obj.minPosition(2) ...
                    || x(3) > obj.maxPosition(3) || x(3) < obj.minPosition(3)
                %warning('Ekf::checkSanity: 3D position out of bounds!');
                obj.stateSanityMask(1) = 1;
            end
            
            % Check 3D attitude
            [roll, pitch, yaw] = quatToRollPitchYaw(x(4:7));
            if roll > obj.maxAttitude(1) || roll < obj.minAttitude(1) ...
                  || pitch > obj.maxAttitude(2) || pitch < obj.minAttitude(2) ...
                  || yaw > obj.maxAttitude(3) || yaw < obj.minAttitude(3)
                  obj.stateSanityMask(2) = 1;
            end
            
            % Check 3D linear velocity
            if x(8) > obj.maxLinVel(1) || x(8) < obj.minLinVel(1) ...
                    || x(9) > obj.maxLinVel(2) || x(9) < obj.minLinVel(2) ...
                    || x(10) > obj.maxLinVel(3) || x(10) < obj.minLinVel(3)
                obj.stateSanityMask(3) = 1;
            end
                
            % Check 3D angular velocity
            if x(11) > obj.maxAngVel(1) || x(11) < obj.minAngVel(1) ...
                    || x(12) > obj.maxAngVel(2) || x(12) < obj.minAngVel(2) ...
                    || x(13) > obj.maxAngVel(3) || x(13) < obj.minAngVel(3)
                obj.stateSanityMask(4) = 1;
            end
            
            % Check 3D linear acceleration
            if x(14) > obj.maxLinAcc(1) || x(14) < obj.minLinAcc(1) ...
                    || x(15) > obj.maxLinAcc(2) || x(15) < obj.minLinAcc(2) ...
                    || x(16) > obj.maxLinAcc(3) || x(16) < obj.minLinAcc(3)
                obj.stateSanityMask(5) = 1;
            end
            
            % Check accelerometer bias
            if x(17) > obj.maxAccBias(1) || x(17) < obj.minAccBias(1) ...
                    || x(18) > obj.maxAccBias(2) || x(18) < obj.minAccBias(2) ...
                    || x(19) > obj.maxAccBias(3) || x(19) < obj.minAccBias(3)
                obj.stateSanityMask(6) = 1;
            end
            
            % Check gyroscope bias
            if x(20) > obj.maxGyrBias(1) || x(20) < obj.minGyrBias(1) ...
                    || x(21) > obj.maxGyrBias(2) || x(21) < obj.minGyrBias(2) ...
                    || x(22) > obj.maxGyrBias(3) || x(22) < obj.minGyrBias(3)
                obj.stateSanityMask(7) = 1;
            end
            
            % --- Check covariance matrix ---
            
            % Reset mask
            obj.covSanityMask = zeros(7, 1);
            
            P = obj.covarianceMatrix;
            
            if P(1, 1) > obj.maxPositionStd(1)^2 ...
                    || P(2, 2) > obj.maxPositionStd(2)^2 ...
                    || P(3, 3) > obj.maxPositionStd(3)^2
                obj.covSanityMask(1) = 1;
            end
            
            if P(4, 4) > obj.maxAttitudeStd(1)^2 ...
                    || P(5, 5) > obj.maxAttitudeStd(2)^2 ...
                    || P(6, 6) > obj.maxAttitudeStd(3)^2
                obj.covSanityMask(2) = 1;
            end
            
            if P(7, 7) > obj.maxLinVelStd(1)^2 ...
                    || P(8, 8) > obj.maxLinVelStd(2)^2 ...
                    || P(9, 9) > obj.maxLinVelStd(3)^2
                obj.covSanityMask(3) = 1;
            end
            
            if P(10, 10) > obj.maxAngVelStd(1)^2 ...
                    || P(11, 11) > obj.maxAngVelStd(2)^2 ...
                    || P(12, 12) > obj.maxAngVelStd(3)^2
                obj.covSanityMask(4) = 1;
            end
            
            if P(13, 13) > obj.maxLinAccStd(1)^2 ...
                    || P(14, 14) > obj.maxLinAccStd(2)^2 ...
                    || P(15, 15) > obj.maxLinAccStd(3)^2
                obj.covSanityMask(5) = 1;
            end
            
            if P(16, 16) > obj.maxAccBiasStd(1)^2 ...
                    || P(17, 17) > obj.maxAccBiasStd(2)^2 ...
                    || P(18, 18) > obj.maxAccBiasStd(3)^2
                obj.covSanityMask(6) = 1;
            end      
            
            if P(19, 19) > obj.maxGyrBiasStd(1)^2 ...
                    || P(20, 20) > obj.maxGyrBiasStd(2)^2 ...
                    || P(21, 21) > obj.maxGyrBiasStd(3)^2
                obj.covSanityMask(7) = 1;
            end 
            
            obj.isSane = (sum(obj.stateSanityMask) + sum(obj.covSanityMask)) == 0;
            
            isGood = obj.isSane;
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
            
            % Check for complex...
            if ~isreal(obj.stateVector)
                %warning('State::checkStateVector: complex detected!');
                for i = 1 : length(obj.stateVector)
                    obj.stateVector(i) = real(obj.stateVector(i));
                end
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
            
            % Check for complex...
            if ~isreal(obj.covarianceMatrix)
                %warning('State::covarianceMatrix: complex detected!');
                for i = 1 : size(obj.covarianceMatrix, 1)
                    for j = 1 : size(obj.covarianceMatrix, 2)
                        obj.covarianceMatrix(i, j) = real(obj.covarianceMatrix(i, j));
                    end
                end
                isGood = false;
            end
        end
    end
end