function [x, P, valid] = stateEstimationFilter(t, z_imu, z_veh, z_gps, z_att, ...
    R_imu, R_veh, R_gps, R_att)
% stateEstimationFilter: function performing one time step of the state estimation filter 
% Args:
% - t (double): current time [s]
% - z_imu (vector [6x1]): [z_acc_x, z_acc_y, z_acc_z, z_gyr_x, z_gyr_y, z_gyr_z] [m/s2] [rad/s]
% - z_veh (vector [7x1]): [w_wheel_fl, w_wheel_fr, w_wheel_rl, w_wheel_rr, delta_f, delta_r, transmission_mode]
% - z_gps (vector): [z_lat, z_lon, z_alt, gps_mode, att_mode, vel_mode, nb_sats]
% - z_att (vector [3x1]): estimated roll, pitch, yaw angles [rad] 
% Return:
% - x (vector [22x1]): estimated state vector
% - P (matrix [21x21]): associated covariance matrix
% - valid (bool): if the estimation is valid

persistent ekf
persistent t_prev
persistent gps_lla_ref
persistent x_prev P_prev

% Initialize global variables
if isempty(ekf)
    ekf = Ekf();
end

if isempty(t_prev)
    t_prev = t;
end

if isempty(x_prev)
    x_prev = zeros(22, 1);
end
if isempty(P_prev)
    P_prev = eye(21, 21);
end

% Conversions
deg2rad = pi / 180.0;
rad2deg = 180.0 / pi;

% disp('t='); disp(t);
% disp('z_imu='); disp(z_imu');
% disp('z_veh='); disp(z_veh');
% disp('z_gps='); disp(z_gps');
% disp('z_att='); disp(z_att');

% Get inputs
wheel_speed_rpm = z_veh(1:4); % [rpm]
steer_angle_rad = z_veh(5:6); % [rad]
transmission_mode = z_veh(7); % 0: sna, 1: park, 2: reverse, 3: neutral, 4: forward

% Compute vehicle speed [m/s]
v_ms = computeVehicleGroundSpeed(wheel_speed_rpm, ekf.wheelRadius);

% Compute vehicle kinematics
[lin_vel, ang_vel] = vehicleKinematics(wheel_speed_rpm, steer_angle_rad, ...
    ekf.wheelRadius, ekf.distFrontAxleToCog, ekf.distRearAxleToCog);

% Check if vehicle is moving
vehicle_at_rest = ((transmission_mode == 1) || (transmission_mode == 3)) ...
    && (abs(v_ms) <= 0.1);

if vehicle_at_rest
    fprintf('Vehicle not moving!\n');
end

% Get and convert GPS data
gps_lla = [z_gps(1) * deg2rad, z_gps(2) * deg2rad, z_gps(3)]; % Convert deg to rad
gps_status = z_gps(4:7); % [pos mode, att mode, vel mode, nb of sat]

% Filter initialization
init_filter = true; % TODO: only init filter when GPS is ready and precision is enough
% init_filter = (gps_status(1) > 0) ...
%     && (gps_status(2) > 0) ...
%     && (gps_status(3) > 0) ...
%     && (gps_status(4) > 4); % TODO

if ~ekf.isInit && init_filter % Filter not initialized AND init conditions are met
    %disp('Initializing the filter...');
    fprintf('Initializing the filter...\n');
    
    % Compute GPS local position in ENU frame
    if isempty(gps_lla_ref)
        gps_lla_ref = gps_lla; % Take first valid GPS measurement as reference for ENU projection
    end
    ecef = lla2ecef(gps_lla); % Convert LLA to ECEF coordinates
    enu = ecef2enu(ecef, gps_lla_ref)'; % Convert ECEF to ENU coordinates
    
    % Initialize state vector
    x = zeros(22, 1);
    x(1:3) = enu; % Init position with GPS [m]
    %x(4:7) = [1.0; 0.0; 0.0; 0.0]; % Init attitude/quaternion
    yaw = z_att(3); % Yaw angle [rad] (from GPS...)
    %fprintf('Init filter with yaw=%3.3f deg\n', yaw * rad2deg);
    x(4:7) = [cos(0.5 * yaw); 0.0; 0.0; sin(0.5 * yaw)]; % Init attitude/quaternion
    x(8:10) = lin_vel; % Init lin vel [m/s] (from vehicle kinematic model)
    x(11:13) = ang_vel; % Init ang vel [rad/s] (from vehicle kinematic model)
    x(14:16) = zeros(3, 1); % Init lin acc [m/s2] (assume no acceleration)
    x(17:19) = zeros(3, 1); % Init acc bias [m/s2] (assume zero bias)
    x(20:22) = zeros(3, 1); % Init gyr bias [rad/s] (assume zero bias)
    
    % Initialize covariance matrix
    P = eye(21, 21);
    % Position uncertainty
    P(1:3, 1:3) = diag([0.3; 0.3; 0.3]).^2; % around 1 m of uncertainty
    % Attitude uncertainty
    P(4:6, 4:6) = diag([0.03; 0.03; 0.06]).^2; % around 5/5/10 deg of uncertainty
    % Linear vel uncertainty
    P(7:9, 7:9) = diag([0.1; 0.1; 0.1]).^2; % around 1 km/h of uncertainty
    % Angular vel uncertainty
    P(10:12, 10:12) = diag([0.01; 0.01; 0.01]).^2;
    % Linear acc uncertainty
    P(13:15, 13:15) = diag([0.01; 0.01; 0.01]).^2;
    % Acc bias uncertainty
    P(16:18, 16:18) = diag([0.15; 0.15; 0.15]).^2;
    % Gyr bias uncertainty
    P(19:21, 19:21) = diag([0.05; 0.05; 0.05]).^2;
    
    % Initialize filter
    ekf.init(t, x, P);
else % Filter is initialized
    % --- Filter prediction ---
    dt = t - t_prev; % EKF period [s]
    assert(dt > 0.0, 'dt=%3.3f must be > 0!', dt);
    %disp('dt='); disp(dt);
    %fprintf('dt=%2.3f\n', dt);
    
    success = ekf.predict(dt);
    if ~success
        %disp('EKF prediction failed!');
        fprintf('EKF prediction failed!\n');
        %warning('EKF prediction failed!');
    end
    % Keep track of prediction
    x_pred = ekf.stateVectorPred;
    P_pred = ekf.covarianceMatrixPred;
    
    % --- Filter correction ---
    
    % Correct with IMU measurements
    z_acc = z_imu(1:3); % Accelero measurement [m/s2]
    z_gyr = z_imu(4:6); % Gyro measurement [rad/s]
    if ~any(isnan(z_imu)) && ~any(isinf(z_imu)) ...
            && (abs(z_acc(1)) < 10.0) && (abs(z_acc(2)) < 10.0) && (abs(z_acc(3)) < 20.0) ...
            && (abs(z_gyr(1)) < 2*pi) && (abs(z_gyr(2)) < 2*pi) && (abs(z_gyr(3)) < 2*pi)
        %disp('IMU correction step');
        %fprintf('IMU correction step: za=(%3.3f, %3.3f, %3.3f) [m/s2], zg=(%3.3f, %3.3f, %3.3f) [rad/s]\n', ...
        %    z_imu(1), z_imu(2), z_imu(3), z_imu(4), z_imu(5), z_imu(6));
        
        maha_thresh_imu = 20.0;%12.59; % chi2(0.05, 6)
        
        ekf.correctImu(z_imu, R_imu, maha_thresh_imu);
    end
    
    % Correct with vehicle measurements
    if vehicle_at_rest % Vehicle not moving
        %disp('Vehicle at rest correction step');
        %fprintf('Vehicle at rest correction step\n');
        % Zero update velocity (ZUV)
        z_zuv = zeros(6, 1);
        R_zuv = diag([0.01; 0.01; 0.01; 0.001; 0.001; 0.001]).^2; % Around 0.1 km/h 3*std, Around 0.2 deg/s 3*std
        maha_thresh_zuv = 25;% 12.59; % chi2(0.05, 6)
        
        ekf.correctLinAngVel(z_zuv, R_zuv, maha_thresh_zuv);
        
        % Zero update acceleration (ZUA)
        z_zua = zeros(3, 1);
        R_zua = diag([0.0033; 0.0033; 0.0033]).^2; % Around 0.01 m/s2 3*std
        maha_thresh_zua = 15;% 7.81; % chi2(0.05, 3)
        
        ekf.correctLinAcc(z_zua, R_zua, maha_thresh_zua);
        
        % Zero update position (ZUP)
        z_zup = x_prev(1:3);
        R_zup = P_prev(1:3, 1:3);
        maha_thresh_zup = 15;% 7.81; % chi2(0.05, 3)
        
        %ekf.correctPosition(z_zup, R_zup, maha_thresh_zup);
        
        % Zero update attitude
        % TODO...
    else % Vehicle moving
        %disp('Vehicle moving correction step');
        %fprintf('Vehicle moving correction step\n');
        
        if lin_vel(2) > 1.0
            disp('vy > 1.0')
            %warning('vy > 1.0');
            %fprintf('vy > 1.0\n');
            lin_vel(2) = 1.0;
        end
        
        z_vel = [lin_vel; ang_vel];

        maha_thresh_vel = 1000.0;% 12.59; % chi2(0.05, 6)
        
        ekf.correctLinAngVel(z_vel, R_veh, maha_thresh_vel);
    end
    
    % Correct with GPS measurements
    if ~isempty(gps_lla_ref) && (gps_status(1) > 0) && (gps_status(4) > 4)
        %disp('GPS correction step');
        %fprintf('GPS correction step\n');
        
        ecef = lla2ecef(gps_lla);
        enu = ecef2enu(ecef, gps_lla_ref)';
        
        maha_thresh_gps = 15;%7.81; % chi2(0.05, 3)

        std_xx = sqrt(R_gps(1, 1));
        std_yy = sqrt(R_gps(2, 2));
        if (3.0 * std_xx < 1.0) && (3.0 * std_yy < 1.0)
            ekf.correctPosition(enu, R_gps, maha_thresh_gps);
        end
    end
    
    % Correct with RT2000 yaw measurement
    %if gps_status(2) > 0 && gps_status(2) < 11
    if true
        %fprintf('Yaw correction step\n');
        
        z_yaw = z_att(3); % Yaw angle [rad]
        
        R_yaw = R_att(3, 3); % Yaw uncertainty

        maha_thresh_yaw = 3.84; % chi2(0.05, 1)
        
        std_yaw = sqrt(R_yaw);

        if (3.0 * std_yaw) < (10.0 * pi / 180.0)
            %fprintf('Correction step with yaw=%3.3f deg\n', z_yaw * rad2deg);
            ekf.correctYaw(z_yaw, R_yaw, maha_thresh_yaw);
        end
    end
    
    % Check filter sanity
    valid = ekf.checkSanity();
    
    if ~valid
        disp('EKF sanity check failed, resetting filter!');
        %fprintf('EKF sanity check failed, resetting filter!\n');
        %warning('EKF sanity check failed, resetting filter!');
        
        ekf.reset();
    end
end

% Set outputs
x = ekf.stateVector;
P = ekf.covarianceMatrix;
valid = ekf.isSane && ekf.isInit;

% Remember previous time [s]
t_prev = t;

% Remember previous state vector and covariance matrix
x_prev = x;
P_prev = P;

end
