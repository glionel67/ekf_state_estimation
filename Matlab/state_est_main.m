% Main script to test the EKF state estimation using Simulink or Matlab
clear all;
close all;
clc;

%% Useful folders
addpath(genpath('D:\Workspace\StateEstimationFilter'));
addpath(genpath('D:\Workspace\StateEstimationFilter\Data'));

%% Options
plot_data = false;
use_simulink = false; % Options: true or false

%% Parameters
deg2rad = pi / 180.0;
rad2deg = 180.0 / pi;
kmh2ms = 1.0 / 3.6;
ms2kmh = 3.6;

% incertitude = 3 * \sigma --> \sigma = incertitude / 3 --> cov = \sigma^2
R_imu = diag([3.0e-2; 3.0e-2; 5.0e-2; 3.0e-3; 3.0e-3; 3.0e-3]).^2; % ax, ay, az, gx, gy, gz
R_veh = diag([0.5/3.0; 0.5/3.0; 0.2/3.0; 0.5*deg2rad/3.0; 0.5*deg2rad/3.0; 0.5*deg2rad/3.0]).^2; % w_fl, w_fr, w_rl, w_rr, d_f, d_r
%R_gps = diag([0.5/3.0; 0.5/3.0; 1.0/3.0]).^2;

%% Load data
data_filename = 'D:\Workspace\StateEstimationFilter\Data\20220912_140003_light.mat';
%data_filename = fullfile(path, filename);

try
    data = load(data_filename);
catch
    error('Failed to load data from file: %s', data_filename);
end

%% Get data

if ~isfield(data, 'Time')
    error('ERROR: no Time signal');
end
if ~isfield(data, 'Time_bis') % For CAN 5 with CSM only
    warning('Creating fake --Time_bis-- variable');
    data.Time_bis = data.Time;
end

nData = length(data.Time);
time_s = data.Time; % [s]
time_min = time_s / 60; % Convert time in [s] to [min]

%% Select range
%rng = 1:nData;
%rng = (time_min > 70) & (time_min < 80); % Vehicle at rest: used to determine IMU noise std
%rng = (time_min > 87) & (time_min < 90); % Vehicle moving: used to test EKF
rng = (time_min > 145) & (time_min < 147);

time_s = data.Time(rng); % [s]
time_min = time_s / 60; % Convert time in [s] to [min]

nData = length(time_s);

%% Extract data

% Motor speeds [rpm]
motor_speed_rpm = [data.MOTOR_SPEEDS_MotorSpeedFL(rng), data.MOTOR_SPEEDS_MotorSpeedFR(rng), ...
    data.MOTOR_SPEEDS_MotorSpeedRL(rng), data.MOTOR_SPEEDS_MotorSpeedRR(rng)]; % [rpm]

% Steering angles [rad]
steer_angle_rad = [data.VEHICLE_DIRECTION_ANGLES_VehicleDirectionAngleFront(rng), ...
    data.VEHICLE_DIRECTION_ANGLES_VehicleDirectionAngleRear(rng)]; % [rad]

vehicle_speed_kmh = data.VEHICLE_STATE_1_VehicleSpeed(rng); % [km/h]
vehicle_mode = data.VEHICLE_STATE_1_VehicleMode(rng);
transmission_mode = data.VEHICLE_STATE_1_TransmissionMode(rng); % 0: sna, 1: park, 2: reverse, 3: neutral, 4: forward

% Accelero measurements
acc_px4_ms2 = [data.IMU_ACC_VehicleLinearAccelX(rng), data.IMU_ACC_VehicleLinearAccelY(rng), data.IMU_ACC_VehicleLinearAccelZ(rng)];
acc_rt_ms2 = [data.AccelVehicle_AccelX(rng), data.AccelVehicle_AccelY(rng), data.AccelVehicle_AccelZ(rng)];
acc_level_rt_ms2 = [data.AccelLevel_AccelForward(rng), data.AccelLevel_AccelLateral(rng), data.AccelLevel_AccelDown(rng)];
% Gyro measurements
gyr_px4_rads = [data.IMU_GYR_VehicleRollRate(rng), data.IMU_GYR_VehiclePitchRate(rng), data.IMU_GYR_VehicleYawRate(rng)];
gyr_rt_rads = [data.RateVehicle_AngRateX(rng), data.RateVehicle_AngRateY(rng), data.RateLevel_AngRateDown(rng)] * deg2rad;
gyr_level_rt_rads = [data.RateLevel_AngRateForward(rng), data.RateLevel_AngRateLateral(rng), data.RateLevel_AngRateDown(rng)] * deg2rad;

%acc_bias_rt_ms2 = [data.AccelBias_AccelBias1(rng), data.AccelBias_AccelBias2(rng), data.AccelBias_AccelBias3(rng)];

vel_level_rt_ms = [data.VelocityLevel_VelForward(rng), data.Velocity_VelEast(rng), data.Velocity_VelDown(rng)];

roll_rt_rad = wrapAngle(data.HeadingPitchRoll_AngleRoll(rng) * deg2rad);
pitch_rt_rad = wrapAngle(data.HeadingPitchRoll_AnglePitch(rng) * deg2rad);
yaw_rt_rad = wrapAngle(data.HeadingPitchRoll_AngleHeading(rng) * deg2rad);
attitude_rt_rad = [roll_rt_rad, pitch_rt_rad, yaw_rt_rad];

gps_lla_px4 = [data.GPS_LAT_LON_Longitude(rng), data.GPS_LAT_LON_Latitude(rng), data.GPS_ALT_FIX_SAT_Altitude(rng)];
gps_lla_rt = [data.LatitudeLongitude_PosLon(rng), data.LatitudeLongitude_PosLat(rng), data.Altitude_PosAlt(rng)];

gps_status_px4 = [data.GPS_ALT_FIX_SAT_FixType(rng), data.GPS_ALT_FIX_SAT_NumberOfSatellites(rng)];

if isfield(data, 'GpsStatus_GpsPosMode')
    gps_status_rt = [data.GpsStatus_GpsPosMode(rng), data.GpsStatus_GpsAttMode(rng), data.GpsStatus_GpsVelMode(rng), data.GpsStatus_GpsNumSats(rng)];
else
    %fprintf('GpsStatus_GpsPosMode message not available!\n');
    warning('GpsStatus_GpsPosMode message not available!');
    gps_status_rt = ones(nData, 4);
end
if isfield(data, 'PosNEDStdev_PosEastStdev')
    pos_cov_rt = [data.PosNEDStdev_PosEastStdev(rng), data.PosNEDStdev_PosNorthStdev(rng), data.PosNEDStdev_PosDownStdev(rng)].^2;
else
    %fprintf('PosNEDStdev_PosEastStdev message not available!\n');
    warning('PosNEDStdev_PosEastStdev message not available!');
    pos_cov_rt = 0.01^2 * ones(nData, 4);
end
if isfield(data, 'AngleStdev_AngleRollStdev')
    att_cov_rt = [data.AngleStdev_AngleRollStdev(rng), data.AngleStdev_AnglePitchStdev(rng), data.AngleStdev_AngleHeadingStdev(rng)].^2;
else
    %fprintf('AngleStdev_AngleRollStdev message not available!\n');
    warning('AngleStdev_AngleRollStdev message not available!');
    att_cov_rt = 0.01^2 * ones(nData, 4);
end

% Plot data (optional)
if plot_data
figure;
subplot(3,1,1);
plot(time_min, gyr_rt_rads(:,1), 'b-');
ylabel('Ang vel X [rad/s]');
title('IMUs gyroscope');
subplot(3,1,2);
plot(time_min, gyr_rt_rads(:,2), 'b-');
ylabel('Ang vel Y [rad/s]');
subplot(3,1,3);
plot(time_min, gyr_rt_rads(:,3), 'b-');
ylabel('Ang vel Z [rad/s]');
xlabel('Time [s]');

figure;
subplot(3,1,1);
plot(time_min, acc_rt_ms2(:,1), 'b-'); hold on;
ylabel('Acc X [m/s2]');
title('IMUs accelerometer');
subplot(3,1,2);
plot(time_min, acc_rt_ms2(:,2), 'b-');
ylabel('Acc Y [m/s2]');
subplot(3,1,3);
plot(time_min, acc_rt_ms2(:,3), 'b-');
ylabel('Acc Z [m/s2]');
xlabel('Time [min]');

figure;
plot(time_min, data.Velocity_Speed2D(rng) * ms2kmh, 'b-'); hold on;
plot(time_min, data.VEHICLE_STATE_1_VehicleSpeed(rng), 'r--');
legend('Speed2D', 'Wheel based');
ylabel('Velocity [km/h]');
xlabel('Time [min]');

figure; hold on;
plot(time_min, motor_speed_rpm(:,1), 'r-');
plot(time_min, motor_speed_rpm(:,2), 'b--');
plot(time_min, motor_speed_rpm(:,3), 'g:');
plot(time_min, motor_speed_rpm(:,4), 'c-.');
legend('FL', 'FR', 'RL', 'RR');
ylabel('Wheel speed [rpm]');
xlabel('Time [min]');
title('Vehicle motor/wheel speeds');

figure; hold on;
plot(time_min, steer_angle_rad(:,1), 'r-');
plot(time_min, steer_angle_rad(:,2), 'b--');
legend('Front', 'Rear');
ylabel('Steering angle [rad]');
xlabel('Time [min]');
title('Vehicle steering angles');

figure;
subplot(3,1,1);
plot(time_min, attitude_rt_rad(:,1) * rad2deg, 'b-');
ylabel('Roll [deg]');
title('Vehicle attitude');
subplot(3,1,2);
plot(time_min, attitude_rt_rad(:,2) * rad2deg, 'b-');
ylabel('Pitch [deg]');
subplot(3,1,3);
plot(time_min, attitude_rt_rad(:,3) * rad2deg, 'b-');
ylabel('Yaw [deg]');
xlabel('Time [min]');

figure;
plot(gps_lla_rt(:,1), gps_lla_rt(:,2), 'r-'); hold on;
plot(gps_lla_px4(:,1), gps_lla_px4(:,2), 'b:');
legend('RT', 'PX4');
xlabel('Longitude [deg]');
ylabel('Latitude [deg]');
title('GPS trajectory');

figure;
plot(time_min, gps_lla_rt(:,3), 'r-'); hold on;
plot(time_min, gps_lla_px4(:,3), 'b:');
legend('RT', 'PX4');
xlabel('Time [min]');
ylabel('Altitude [m]');
title('GPS altitude');

figure; hold on;
plot(time_min, gps_status_rt(:,1), 'ro');
plot(time_min, gps_status_rt(:,2), 'g.');
plot(time_min, gps_status_rt(:,3), 'b+');
legend('RT pos', 'RT att', 'RT vel');
xlabel('Time [min]');
ylabel('Status');
title('GPS status');

figure;
subplot(3,1,1);
plot(time_min, sqrt(pos_cov_rt(:,1)), 'b-');
ylabel('Pos X [m]');
title('Pos std');
subplot(3,1,2);
plot(time_min, sqrt(pos_cov_rt(:,2)), 'b-');
ylabel('Pos Y [m]');
subplot(3,1,3);
plot(time_min, sqrt(pos_cov_rt(:,3)), 'b-');
ylabel('Pos Z [m]');
xlabel('Time [s]');

figure;
subplot(3,1,1);
plot(time_min, sqrt(att_cov_rt(:,1)), 'b-');
ylabel('Att X [rad]');
title('Att std');
subplot(3,1,2);
plot(time_min, sqrt(att_cov_rt(:,2)), 'b-');
ylabel('Att Y [rad]');
subplot(3,1,3);
plot(time_min, sqrt(att_cov_rt(:,3)), 'b-');
ylabel('Att Z [rad]');
xlabel('Time [s]');

end

%% Simulation

if use_simulink
    %% Open simulink model
    simulinkModelName = 'StateEstimationFilter'; %'StateEstimationFilter.slx';
    open_system(simulinkModelName);
    
    %% Run simulation
    diary
    simInModel = Simulink.SimulationInput(simulinkModelName);
    simInModel = setModelParameter(simInModel, 'StartTime', num2str(time_s(1)), 'StopTime', num2str(time_s(end)));
    
    simOutModel = sim(simInModel);
    stateVector = simOutModel.stateVector;
    covarianceMatrix = simOutModel.covarianceMatrix;
    isValid = simOutModel.isValid;
    diary OFF
    
    %% Extract results from simulation
    time_res_s = time_s(1) + stateVector.Time;
    nResData = length(time_s);
    position = stateVector.Data(:,1:3);
    quaternion = stateVector.Data(:,4:7);
    attitude = zeros(nResData, 3);
    for n = 1 : nResData
        [roll, pitch, yaw] = quatToRollPitchYaw(quaternion(n, 1:4));
        attitude(n, 1:3) = [roll, pitch, yaw];
    end
    linVel = stateVector.Data(:,8:10);
    angVel = stateVector.Data(:,11:13);
    linAcc = stateVector.Data(:,14:16);
    accBias = stateVector.Data(:,17:19);
    gyrBias = stateVector.Data(:,20:22);
    
    posCov = permute(covarianceMatrix.Data(1:3, 1:3, :), [3, 1, 2]);
    attCov = permute(covarianceMatrix.Data(4:6, 4:6, :), [3, 1, 2]);
    linVelCov = permute(covarianceMatrix.Data(7:9, 7:9, :), [3, 1, 2]);
    angVelCov = permute(covarianceMatrix.Data(10:12, 10:12, :), [3, 1, 2]);
    linAccCov = permute(covarianceMatrix.Data(13:15, 13:15, :), [3, 1, 2]);
    accBiasCov = permute(covarianceMatrix.Data(16:18, 16:18, :), [3, 1, 2]);
    gyrBiasCov = permute(covarianceMatrix.Data(19:21, 19:21, :), [3, 1, 2]);
    
    valid = isValid.data;
    
else
    diary
    nb_iters = length(time_s);
    nb_states = 22;
	stateVector = zeros(nb_iters, nb_states);
    covarianceMatrix = cell(nb_iters, 1);
    isValid = zeros(nb_iters, 1);
    
    lin_vel_veh = zeros(nb_iters, 3);
    ang_vel_veh = zeros(nb_iters, 3);
    pos_veh = zeros(nb_iters, 3);
    att_veh = zeros(nb_iters, 3);
    
    for n = 1 : nb_iters
        t = time_s(n);
        z_imu = [acc_rt_ms2(n, :)'; gyr_rt_rads(n, :)'];
        z_veh = [motor_speed_rpm(n, :)'; steer_angle_rad(n, :)'; transmission_mode(n)];
        z_gps = [gps_lla_rt(n, :)'; gps_status_rt(n, :)'];
        z_att = attitude_rt_rad(n, :)';
        
        %R_imu = ...; % Fixed value
        %R_veh = ...; % Fixed value
        R_gps = diag(pos_cov_rt(n, :));
        R_att = diag(att_cov_rt(n, :));
        
        % Perform one step of Kalman filter
        [x, P, valid] = stateEstimationFilter(t, z_imu, z_veh, z_gps, z_att, ...
            R_imu, R_veh, R_gps, R_att);
        
        stateVector(n, :) = x';
        covarianceMatrix{n, 1} = P;
        isValid(n, 1) = valid;
        
        % Vehicle odometry only
        [v, w] = vehicleKinematics(motor_speed_rpm(n, :), steer_angle_rad(n, :), 0.3575, 1.115, 1.135);
        lin_vel_veh(n,:) = v;
        ang_vel_veh(n,:) = w;
        if n == 1
            pos_veh(n,:) = stateVector(n, 1:3);
            [roll, pitch, yaw] = quatToRollPitchYaw(stateVector(n, 4:7)');
            att_veh(n,:) = [roll, pitch, yaw];
        else
            dt = time_s(n) - time_s(n - 1);
            R = rollPitchYawToRotMat(att_veh(n - 1,1), att_veh(n - 1,2), att_veh(n - 1,3));
            pos_veh(n,:) = pos_veh(n - 1,:) + (R * lin_vel_veh(n, :)' * dt)';
            att_veh(n,:) = wrapAngle(att_veh(n - 1,:) + ang_vel_veh(n, :) * dt);
        end
        
    end
    diary off
    
    %% Extract results from simulation
    time_res_s = time_s;
    position = stateVector(:,1:3);
    quaternion = stateVector(:,4:7);
    attitude = zeros(nb_iters, 3);
    linVel = stateVector(:,8:10);
    angVel = stateVector(:,11:13);
    linAcc = stateVector(:,14:16);
    accBias = stateVector(:,17:19);
    gyrBias = stateVector(:,20:22);
    
    posCov = zeros(nb_iters, 3, 3);
    attCov = zeros(nb_iters, 3, 3);
    linVelCov = zeros(nb_iters, 3, 3);
    angVelCov = zeros(nb_iters, 3, 3);
    linAccCov = zeros(nb_iters, 3, 3);
    accBiasCov = zeros(nb_iters, 3, 3);
    gyrBiasCov = zeros(nb_iters, 3, 3);
    
    valid = isValid;
    
    for n = 1 : nb_iters
        [roll, pitch, yaw] = quatToRollPitchYaw(quaternion(n, 1:4));
        attitude(n, 1:3) = [roll, pitch, yaw];
        
        posCov(n,:,:) = covarianceMatrix{n, 1}(1:3, 1:3);
        attCov(n,:,:) = covarianceMatrix{n, 1}(4:6, 4:6);
        linVelCov(n,:,:) = covarianceMatrix{n, 1}(7:9, 7:9);
        angVelCov(n,:,:) = covarianceMatrix{n, 1}(10:12, 10:12);
        linAccCov(n,:,:) = covarianceMatrix{n, 1}(13:15, 13:15);
        accBiasCov(n,:,:) = covarianceMatrix{n, 1}(16:18, 16:18);
        gyrBiasCov(n,:,:) = covarianceMatrix{n, 1}(19:21, 19:21);
    end
end

nb_of_reset = sum(valid == 0);
fprintf('Number of filter reset: %d\n', nb_of_reset);

%% Plot result
figure;
plot(time_min, valid, 'b.');
xlabel('Time [min]');
ylabel('Valid [0/1]');
title('Estimation validity');

% Position
lla = [gps_lla_rt(:,1) * deg2rad, gps_lla_rt(:,2) * deg2rad, gps_lla_rt(:,3)];
ecef = lla2ecef(lla);
enu = ecef2enu(ecef, lla(1,:));
figure;
plot(enu(:,1), enu(:,2), 'r-'); hold on;
plot(position(:,1), position(:,2), 'b--'); hold on;
plot(position(1,1), position(1,2), 'ks');
plot(position(end,1), position(end,2), 'gd');
plot(pos_veh(:,1), pos_veh(:,2), 'c:');
r = 1:100:length(position);
quiver(position(r,1), position(r,2), 2*cos(attitude(r,3)), 2*sin(attitude(r,3)), 0, 'k-');
quiver(enu(r,1), enu(r,2), 2*cos(attitude_rt_rad(r,3)), 2*sin(attitude_rt_rad(r,3)), 0, 'k-');
legend('RT2000', 'EKF');
xlabel('X [m]');
ylabel('Y [m]');
title('Position');

% Attitude
figure;
subplot(3,1,1);
plot(time_min, attitude_rt_rad(:,1) * rad2deg, 'r-'); hold on;
plot(time_min, attitude(:,1) * rad2deg, 'b--');
plot(time_min, att_veh(:,1) * rad2deg, 'c:');
legend('RT2000', 'EKF');
ylabel('Roll [deg]');
title('Attitude');
subplot(3,1,2);
plot(time_min, attitude_rt_rad(:,2) * rad2deg, 'r-'); hold on;
plot(time_min, attitude(:,2) * rad2deg, 'b--');
plot(time_min, att_veh(:,2) * rad2deg, 'c:');
legend('RT2000', 'EKF');
ylabel('Pitch [deg]');
subplot(3,1,3);
plot(time_min, attitude_rt_rad(:,3) * rad2deg, 'r-');  hold on;
plot(time_min, attitude(:,3) * rad2deg, 'b--');
plot(time_min, att_veh(:,3) * rad2deg, 'c:');
legend('RT2000', 'EKF');
ylabel('Yaw [deg]');
xlabel('Time [min]');

% Linear velocity
figure;
subplot(3,1,1);
plot(time_min, data.VelocityLevel_VelForward(rng), 'r-'); hold on;
plot(time_min, linVel(:,1), 'b:');
plot(time_min, lin_vel_veh(:,1), 'c:');
legend('RT2000', 'EKF');
ylabel('Vel x [m/s]');
subplot(3,1,2);
plot(time_min, data.VelocityLevel_VelLateral(rng), 'r-'); hold on;
plot(time_min, linVel(:,2), 'b:');
plot(time_min, lin_vel_veh(:,2), 'c:');
legend('RT2000', 'EKF');
ylabel('Vel y [m/s]');
subplot(3,1,3);
plot(time_min, data.Velocity_VelDown(rng), 'r-'); hold on;
plot(time_min, linVel(:,3), 'b:');
plot(time_min, lin_vel_veh(:,3), 'c:');
legend('RT2000', 'EKF');
ylabel('Vel z [m/s]');
xlabel('Time [s]');

% Angular velocity
figure;
subplot(3,1,1);
plot(time_min, data.RateLevel_AngRateForward(rng), 'r-'); hold on;
plot(time_min, angVel(:,1) * rad2deg, 'b--');
plot(time_min, ang_vel_veh(:,1) * rad2deg, 'c:');
legend('RT2000', 'EKF');
ylabel('Vel X [deg/s]');
title('Angular velocity');
subplot(3,1,2);
plot(time_min, data.RateLevel_AngRateLateral(rng), 'r-'); hold on;
plot(time_min, angVel(:,2) * rad2deg, 'b--');
plot(time_min, ang_vel_veh(:,2) * rad2deg, 'c:');
legend('RT2000', 'EKF');
ylabel('Vel Y [deg/s]');
subplot(3,1,3);
plot(time_min, data.RateLevel_AngRateDown(rng), 'r-');  hold on;
plot(time_min, angVel(:,3) * rad2deg, 'b--');
plot(time_min, ang_vel_veh(:,3) * rad2deg, 'c:');
legend('RT2000', 'EKF');
ylabel('Vel Z [deg/s]');
xlabel('Time [min]');

% Linear acceleration
figure;
subplot(3,1,1);
plot(time_min, data.AccelLevel_AccelForward(rng), 'r-'); hold on;
plot(time_min, linAcc(:,1), 'b--');
legend('RT2000', 'EKF');
ylabel('Acc X [m/s2]');
title('Linear acceleration');
subplot(3,1,2);
plot(time_min, data.AccelLevel_AccelLateral(rng), 'r-'); hold on;
plot(time_min, linAcc(:,2), 'b--');
legend('RT2000', 'EKF');
ylabel('Acc Y [m/s2]');
subplot(3,1,3);
plot(time_min, data.AccelLevel_AccelDown(rng), 'r-');  hold on;
plot(time_min, linAcc(:,3), 'b--');
legend('RT2000', 'EKF');
ylabel('Acc Z [m/s2]');
xlabel('Time [min]');

% Acc bias
figure;
subplot(3,1,1);
plot(time_min, accBias(:,1), 'b--');
legend('EKF');
ylabel('X [m/s2]');
title('Accelero bias');
subplot(3,1,2);
plot(time_min, accBias(:,2), 'b--');
legend('EKF');
ylabel('Y [m/s2]');
subplot(3,1,3);
plot(time_min, accBias(:,3), 'b--');
legend('EKF');
ylabel('Z [m/s2]');
xlabel('Time [min]');

% Gyr bias
figure;
subplot(3,1,1);
plot(time_min, gyrBias(:,1), 'b--');
legend('EKF');
ylabel('X [rad/s]');
title('Gyro bias');
subplot(3,1,2);
plot(time_min, gyrBias(:,2), 'b--');
legend('EKF');
ylabel('Y [rad/s]');
subplot(3,1,3);
plot(time_min, gyrBias(:,3), 'b--');
legend('EKF');
ylabel('Z [rad/s]');
xlabel('Time [min]');

% Position Covariance
figure; hold on;
plot(time_min, sqrt(posCov(:,1,1)), 'r-');
plot(time_min, sqrt(posCov(:,2,2)), 'b--');
plot(time_min, sqrt(posCov(:,3,3)), 'g:');
legend('x', 'y', 'z');
xlabel('Time [min]');
ylabel('Pos std [m]');
title('Position filter std');

% Attitude covariance
figure; hold on;
plot(time_min, sqrt(attCov(:,1,1)), 'r-');
plot(time_min, sqrt(attCov(:,2,2)), 'b--');
plot(time_min, sqrt(attCov(:,3,3)), 'g:');
legend('x', 'y', 'z');
xlabel('Time [min]');
ylabel('Att std [rad]');
title('Attitude filter std');

% Lin vel covariance
figure; hold on;
plot(time_min, sqrt(linVelCov(:,1,1)), 'r-');
plot(time_min, sqrt(linVelCov(:,2,2)), 'b--');
plot(time_min, sqrt(linVelCov(:,3,3)), 'g:');
legend('x', 'y', 'z');
xlabel('Time [min]');
ylabel('Lin vel std [m/s]');
title('Linear velocity filter std');

% Ang vel covariance
figure; hold on;
plot(time_min, sqrt(angVelCov(:,1,1)), 'r-');
plot(time_min, sqrt(angVelCov(:,2,2)), 'b--');
plot(time_min, sqrt(angVelCov(:,3,3)), 'g:');
legend('x', 'y', 'z');
xlabel('Time [min]');
ylabel('Ang vel std [rad/s]');
title('Angular velocity filter std');

% Lin acc covariance
figure; hold on;
plot(time_min, sqrt(linAccCov(:,1,1)), 'r-');
plot(time_min, sqrt(linAccCov(:,2,2)), 'b--');
plot(time_min, sqrt(linAccCov(:,3,3)), 'g:');
legend('x', 'y', 'z');
xlabel('Time [min]');
ylabel('Lin acc std [m/s2]');
title('Linear acceleration filter std');

% Acc bias covariance
figure; hold on;
plot(time_min, sqrt(accBiasCov(:,1,1)), 'r-');
plot(time_min, sqrt(accBiasCov(:,2,2)), 'b--');
plot(time_min, sqrt(accBiasCov(:,3,3)), 'g:');
legend('x', 'y', 'z');
xlabel('Time [min]');
ylabel('Acc bias std [m/s2]');
title('Accelero bias filter std');

% Gyr bias covariance
figure; hold on;
plot(time_min, sqrt(gyrBiasCov(:,1,1)), 'r-');
plot(time_min, sqrt(gyrBiasCov(:,2,2)), 'b--');
plot(time_min, sqrt(gyrBiasCov(:,3,3)), 'g:');
legend('x', 'y', 'z');
xlabel('Time [min]');
ylabel('Gyr bias std [rad/s]');
title('Gyro bias filter std');