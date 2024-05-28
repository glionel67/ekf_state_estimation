% Main script to test the EKF state estimation using Simulink or Matlab
clear all;
close all;
clc;

%% Useful folders


%% Options
use_simulink = false; % Options: true or false
plot_data = false;

%% Conversions
deg2rad = pi / 180.0;
rad2deg = 180.0 / pi;
rads2rpm = 30.0 / pi;
rpm2rads = pi / 30.0;
kmh2ms = 1.0 / 3.6;
ms2kmh = 3.6;

%% Parameters
R_imu = diag([3.0e-2; 3.0e-2; 5.0e-2; 3.0e-3; 3.0e-3; 3.0e-3]).^2; % ax, ay, az, gx, gy, gz
R_veh = diag([1.0; 1.0; 1.0; 0.5*deg2rad/3.0; 0.5*deg2rad/3.0; 0.5*deg2rad/3.0]).^2; % w_fl, w_fr, w_rl, w_rr, d_f, d_r
R_gps = diag([0.1; 0.1; 1.0]).^2;

% Encoder
encoder_resolution = 4096.0; % Encoder resolution
encoder_left_wheel_diameter = 0.623479; % [m]
encoder_right_wheel_diameter = 0.622806; % [m]
encoder_wheel_base = 1.52439; % [m]

%% Choose sequence
seq_dir = 'D:\Datasets\urban27-dongtan';

%% Import data

% Global pose [timestamp, P(0,0), P(0,1), P(0,2), P(0,3), P(1,0), P(1,1), P(1,2), P(1,3), P(2,0), P(2,1), P(2,2), P(2,3)]
% P is 4X4 transformation matrix at each time. Rotation matrix and translation matrix is stored in vector form.
filename = fullfile(seq_dir, 'global_pose.csv');
global_pose = readmatrix(filename, 'Delimiter', ',');
%global_pose_ts = readmatrix(filename, 'Delimiter', ',', 'OutputType', 'uint64', 'Range', 'A:A');

%global_pose_xyz = [global_pose(:,5), global_pose(:,9), global_pose(:,13)];

dt = diff(global_pose(:,1)) * 1e-9;
dt_mean = mean(dt);
dt_median = median(dt);
fprintf('Global pose dt mean=%3.4f s, median=%3.4f s\n', dt_mean, dt_median);

figure;
plot(global_pose(:,5), global_pose(:, 9), 'b-');
xlabel('X [m]'); ylabel('Y [m]'); grid on;
title('Baseline trajectory estimation');

% Encoder [timestamp, left count, right count]
filename = fullfile(seq_dir, 'sensor_data/encoder.csv');
encoder = readmatrix(filename, 'Delimiter', ',');
% Compute wheel speeds
enc_diff = [zeros(1, 3); diff(encoder)];
enc_diff(:,1) = enc_diff(:,1) * 1e-9; % Convert to [s]
enc_diff(:,2:3) = 2 * pi * enc_diff(:,2:3) / encoder_resolution; %Convert to [rad]
enc_speed = [encoder(:,1), ...
    0.5 * encoder_left_wheel_diameter * enc_diff(:,2) ./ enc_diff(:,1), ...
    0.5 * encoder_right_wheel_diameter * enc_diff(:,3) ./ enc_diff(:,1)];
enc_speed(1,2:3) = zeros(1, 2);

dt = enc_diff(:,1);
dt_mean = mean(dt);
dt_median = median(dt);
fprintf('Encoder dt mean=%3.4f s, median=%3.4f s\n', dt_mean, dt_median);

% GPS: [timestamp, latitude, longitude, altitude, 9-tuple vector (position covariance)]
filename = fullfile(seq_dir, 'sensor_data/gps.csv');
gps = readmatrix(filename, 'Delimiter', ',');
%lla = gps(:,2:4);
%lla_std = [gps(:,5), gps(:,9), gps(:,13)];

dt = diff(gps(:,1)) * 1e-9;
dt_mean = mean(dt);
dt_median = median(dt);
fprintf('GPS dt mean=%3.4f s, median=%3.4f s\n', dt_mean, dt_median);

% VRS GPS
% Ver1: [timestamp, lat, long, x_UTM, y_UTM, alt, fix, nb of sat, 
% horiz precision, lat std, long std, alt std, heading valid flag, 
% mag global heading, speed in knot, speed in km, GNVTG mode]
% Ver2: Ver1 + ortometric altitude]
% Fix: 4=fix, 5=float, 1=normal
filename = fullfile(seq_dir, 'sensor_data/vrs_gps.csv');
vrs_gps = readmatrix(filename, 'Delimiter', ',');
% lla_vrs = [vrs_gps(:, 2:3), vrs_gps(:, 6)];
% lla_vrs_std = [vrs_gps(:, 10:12)];
speed_kmh = vrs_gps(:, 16);

dt = diff(vrs_gps(:,1)) * 1e-9;
dt_mean = mean(dt);
dt_median = median(dt);
fprintf('VRS GPS dt mean=%3.4f s, median=%3.4f s\n', dt_mean, dt_median);

figure;
plot(vrs_gps(:,3), vrs_gps(:,2), 'r-'); hold on;
plot(gps(:,3), gps(:,2), 'b--');
legend('VRS GPS', 'GPS');
xlabel('Lon [deg]'); ylabel('Lat [deg]'); grid on;
title('GPS');

figure;
plot(vrs_gps(:,1)*1e-9, speed_kmh, 'r-'); hold on;
plot(enc_speed(:,1)*1e-9, enc_speed(:,2)*3.6, 'b--');
plot(enc_speed(:,1)*1e-9, enc_speed(:,3)*3.6, 'g:');
xlabel('Time [s]'); ylabel('Speed [km/h]'); grid on;
legend('VRS GPS', 'ENC left', 'ENC right');
title('Vehicle speed comparison');

% IMU
% Ver1: [timestamp, qx, qy, qz, qw, Euler x, Euler y, Euler z]
% Ver2: Ver1 + [gx, gy, gz, ax, ay, az, mx, my, mz]
filename = fullfile(seq_dir, 'sensor_data/xsens_imu.csv');
imu = readmatrix(filename, 'Delimiter', ',');
% quat = imu(:,2:5);
% att = imu(:,6:8);
% gyr = imu(:,9:11);
% acc = imu(:,12:14);
% mag = imu(:,15:17);

dt = diff(imu(:,1)) * 1e-9;
dt_mean = mean(dt);
dt_median = median(dt);
fprintf('IMU dt mean=%3.4f s, median=%3.4f s\n', dt_mean, dt_median);

% FOG (Fiber optic Gyro) [timestamp, delta roll, delta pitch, delta yaw]
filename = fullfile(seq_dir, 'sensor_data/fog.csv');
fog = readmatrix(filename, 'Delimiter', ',');

dt = diff(fog(:,1)) * 1e-9;
dt_mean = mean(dt);
dt_median = median(dt);
fprintf('FOG dt mean=%3.4f s, median=%3.4f s\n', dt_mean, dt_median);

% Altimeter [timestamp, altitude]
filename = fullfile(seq_dir, 'sensor_data/altimeter.csv');
altimeter = readmatrix(filename, 'Delimiter', ',');

dt = diff(altimeter(:,1)) * 1e-9;
dt_mean = mean(dt);
dt_median = median(dt);
fprintf('Altimeter dt mean=%3.4f s, median=%3.4f s\n', dt_mean, dt_median);

% Data [timestamp, data name]
filename = fullfile(seq_dir, 'sensor_data/data_stamp.csv');
data_stamp = readmatrix(filename, 'Delimiter', ',');
t0 = data_stamp(1,1);


%% Resample data
time = imu(:,1);
time_s = time * 1e-9;
time_min = time_s / 60.0;

nData = length(time);

global_pose2 = zeros(length(time), size(global_pose, 2));
for i = 2 : size(global_pose, 2)
    global_pose2(:,i) = interp1(global_pose(:,1), global_pose(:,i), time, 'nearest', 'extrap');
end
global_pose2(:,1) = time;

enc_speed2 = zeros(length(time), size(enc_speed, 2));
for i = 2 : size(enc_speed, 2)
    enc_speed2(:,i) = interp1(enc_speed(:,1), enc_speed(:,i), time, 'nearest', 'extrap');
end
enc_speed2(:,1) = time;

gps2 = zeros(length(time), size(gps, 2));
for i = 2 : size(gps, 2)
    gps2(:,i) = interp1(gps(:,1), gps(:,i), time, 'nearest', 'extrap');
end
gps2(:,1) = time;

vrs_gps2 = zeros(length(time), size(vrs_gps, 2));
for i = 2 : size(vrs_gps, 2)
    vrs_gps2(:,i) = interp1(vrs_gps(:,1), vrs_gps(:,i), time, 'nearest', 'extrap');
end
vrs_gps2(:,1) = time;

fog2 = zeros(length(time), size(fog, 2));
for i = 2 : size(fog, 2)
    fog2(:,i) = interp1(fog(:,1), fog(:,i), time, 'nearest', 'extrap');
end
fog2(:,1) = time;

altimeter2 = zeros(length(time), size(altimeter, 2));
for i = 2 : size(altimeter, 2)
    altimeter2(:,i) = interp1(altimeter(:,1), altimeter(:,i), time, 'nearest', 'extrap');
end
altimeter2(:,1) = time;

%% Select range
rng = 1:nData;
%rng = (time_min > 70) & (time_min < 80); % Vehicle at rest: used to determine IMU noise std

time_s = time_s(rng); % [s]
time_min = time_s / 60; % Convert time in [s] to [min]

nData = length(time_s);

%% Extract data

% Motor speeds [rpm]
motor_speed_rpm = [zeros(nData, 1), zeros(nData, 1), ...
    enc_speed2(rng, 2) * rads2rpm * 2 / encoder_left_wheel_diameter, ...
    enc_speed2(rng, 3) * rads2rpm * 2 / encoder_right_wheel_diameter]; % [rpm]

% Steering angles [rad]
steer_angle_rad = [zeros(nData, 1), ...
    zeros(nData, 1)]; % [rad]

% Vehicle transmission
transmission_mode = 4 * ones(nData, 1); % 0: sna, 1: park, 2: reverse, 3: neutral, 4: forward
mask = (enc_speed2(rng, 2) == 0) & (enc_speed2(rng, 3) == 0);
transmission_mode(mask) = 1;

% Accelero measurements
acc_ms2 = imu(rng, 12:14);

% Gyro measurements
gyr_rads = imu(rng, 9:11);

% Attitude measurements
attitude_rad = zeros(nData, 3);
for n = 1 : nData
    q = [imu(n, 5), imu(2:4)];
    [attitude_rad(n, 1), attitude_rad(n, 2), attitude_rad(n, 3)] = ...
        quatToRollPitchYaw(q);
end

% GPS measurements
gps_lla = gps2(rng, 2:4);
gps_status = [zeros(nData, 1), zeros(nData, 1), zeros(nData, 1), zeros(nData, 1)];
gps_pos_cov = [gps2(rng, 5), gps2(rng, 9), gps2(rng, 13)];
att_cov = zeros(nData, 3);

vrs_gps_lla = [vrs_gps2(:, 2:3), vrs_gps2(:, 6)];
vrs_gps_pos_cov = vrs_gps2(:, 10:12);

% Plot data (optional)
if plot_data
figure;
subplot(3,1,1);
plot(time_min, gyr_rads(:,1), 'b-');
ylabel('Ang vel X [rad/s]');
title('IMUs gyroscope');
subplot(3,1,2);
plot(time_min, gyr_rads(:,2), 'b-');
ylabel('Ang vel Y [rad/s]');
subplot(3,1,3);
plot(time_min, gyr_rads(:,3), 'b-');
ylabel('Ang vel Z [rad/s]');
xlabel('Time [s]');

figure;
subplot(3,1,1);
plot(time_min, acc_ms2(:,1), 'b-'); hold on;
ylabel('Acc X [m/s2]');
title('IMUs accelerometer');
subplot(3,1,2);
plot(time_min, acc_ms2(:,2), 'b-');
ylabel('Acc Y [m/s2]');
subplot(3,1,3);
plot(time_min, acc_ms2(:,3), 'b-');
ylabel('Acc Z [m/s2]');
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
plot(time_min, attitude_rad(:,1) * rad2deg, 'b-');
ylabel('Roll [deg]');
title('Vehicle attitude');
subplot(3,1,2);
plot(time_min, attitude_rad(:,2) * rad2deg, 'b-');
ylabel('Pitch [deg]');
subplot(3,1,3);
plot(time_min, attitude_rad(:,3) * rad2deg, 'b-');
ylabel('Yaw [deg]');
xlabel('Time [min]');

figure;
plot(vrs_gps_lla(:,2), vrs_gps_lla(:,1), 'r-'); hold on;
plot(gps_lla(:,2), gps_lla(:,1), 'b--');
legend('VRS GPS', 'GPS');
xlabel('Longitude [deg]');
ylabel('Latitude [deg]');
title('GPS trajectory');

figure;
plot(time_min, vrs_gps_lla(:,3), 'r-'); hold on;
plot(time_min, gps_lla(:,3), 'b--');
legend('VRS GPS', 'GPS');
xlabel('Time [min]');
ylabel('Altitude [m]');
title('GPS altitude');

figure; hold on;
plot(time_min, gps_status(:,1), 'ro');
plot(time_min, gps_status(:,2), 'g.');
plot(time_min, gps_status(:,3), 'b+');
legend('RT pos', 'RT att', 'RT vel');
xlabel('Time [min]');
ylabel('Status');
title('GPS status');

figure;
subplot(3,1,1); hold on;
plot(time_min, sqrt(vrs_gps_pos_cov(:,1)), 'r-');
plot(time_min, sqrt(gps_pos_cov(:,1)), 'b--');
ylabel('Pos X [m]');
title('Pos std');
subplot(3,1,2); hold on;
plot(time_min, sqrt(vrs_gps_pos_cov(:,2)), 'r-');
plot(time_min, sqrt(gps_pos_cov(:,2)), 'b--');
ylabel('Pos Y [m]');
subplot(3,1,3); hold on;
plot(time_min, sqrt(vrs_gps_pos_cov(:,3)), 'r-');
plot(time_min, sqrt(gps_pos_cov(:,3)), 'b--');
ylabel('Pos Z [m]');
xlabel('Time [s]');

figure;
subplot(3,1,1);
plot(time_min, sqrt(att_cov(:,1)), 'b-');
ylabel('Att X [rad]');
title('Att std');
subplot(3,1,2);
plot(time_min, sqrt(att_cov(:,2)), 'b-');
ylabel('Att Y [rad]');
subplot(3,1,3);
plot(time_min, sqrt(att_cov(:,3)), 'b-');
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
    for n = 1 : nb_iters
        t = time_s(n);
        z_imu = [acc_ms2(n, :)'; gyr_rads(n, :)'];
        z_veh = [motor_speed_rpm(n, :)'; steer_angle_rad(n, :)'; transmission_mode(n)];
        z_gps = [gps_lla(n, :)'; gps_status(n, :)'];
        z_att = attitude_rad(n, :)';
        
        %R_imu = ...; % Fixed value
        %R_veh = ...; % Fixed value
        R_gps = diag(gps_pos_cov(n, :));
        R_att = diag(att_cov(n, :));
        
        [x, P, valid] = stateEstimationFilter(t, z_imu, z_veh, z_gps, z_att, ...
            R_imu, R_veh, R_gps, R_att);
        
        stateVector(n, :) = x';
        covarianceMatrix{n, 1} = P;
        isValid(n, 1) = valid;
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

%% Plot result
figure;
plot(time_min, valid, 'b.');
xlabel('Time [min]');
ylabel('Valid [0/1]');
title('Estimation validity');

% Position
lla = [vrs_gps_lla(:,1) * deg2rad, vrs_gps_lla(:,2) * deg2rad, vrs_gps_lla(:,3)];
ecef = lla2ecef(lla);
enu = ecef2enu(ecef, lla(1,:));
figure;
plot(enu(:,1), enu(:,2), 'r-'); hold on;
plot(position(:,1), position(:,2), 'b--'); hold on;
plot(position(1,1), position(1,2), 'ks');
plot(position(end,1), position(end,2), 'gd');
r = 1:100:length(position);
quiver(position(r,1), position(r,2), 2*cos(attitude(r,3)), 2*sin(attitude(r,3)), 0, 'k-');
quiver(enu(r,1), enu(r,2), 2*cos(attitude_rad(r,3)), 2*sin(attitude_rad(r,3)), 0, 'k-');
legend('RT2000', 'EKF');
xlabel('X [m]');
ylabel('Y [m]');
title('Position');

% Attitude
figure;
subplot(3,1,1);
plot(time_min, attitude_rad(:,1) * rad2deg, 'r-'); hold on;
plot(time_min, attitude(:,1) * rad2deg, 'b--');
legend('RT2000', 'EKF');
ylabel('Roll [deg]');
title('Attitude');
subplot(3,1,2);
plot(time_min, attitude_rad(:,2) * rad2deg, 'r-'); hold on;
plot(time_min, attitude(:,2) * rad2deg, 'b--');
legend('RT2000', 'EKF');
ylabel('Pitch [deg]');
subplot(3,1,3);
plot(time_min, attitude_rad(:,3) * rad2deg, 'r-');  hold on;
plot(time_min, attitude(:,3) * rad2deg, 'b--');
legend('RT2000', 'EKF');
ylabel('Yaw [deg]');
xlabel('Time [min]');

% Linear velocity
figure;
subplot(3,1,1);
plot(time_min, data.VelocityLevel_VelForward(rng), 'r-'); hold on;
plot(time_min, linVel(:,1), 'b:');
legend('RT2000', 'EKF');
ylabel('Vel x [m/s]');
subplot(3,1,2);
plot(time_min, data.VelocityLevel_VelLateral(rng), 'r-'); hold on;
plot(time_min, linVel(:,2), 'b:');
legend('RT2000', 'EKF');
ylabel('Vel y [m/s]');
subplot(3,1,3);
plot(time_min, data.Velocity_VelDown(rng), 'r-'); hold on;
plot(time_min, linVel(:,3), 'b:');
legend('RT2000', 'EKF');
ylabel('Vel z [m/s]');
xlabel('Time [s]');

% Angular velocity
figure;
subplot(3,1,1);
plot(time_min, data.RateLevel_AngRateForward(rng), 'r-'); hold on;
plot(time_min, angVel(:,1) * rad2deg, 'b--');
legend('RT2000', 'EKF');
ylabel('Vel X [deg/s]');
title('Angular velocity');
subplot(3,1,2);
plot(time_min, data.RateLevel_AngRateLateral(rng), 'r-'); hold on;
plot(time_min, angVel(:,2) * rad2deg, 'b--');
legend('RT2000', 'EKF');
ylabel('Vel Y [deg/s]');
subplot(3,1,3);
plot(time_min, data.RateLevel_AngRateDown(rng), 'r-');  hold on;
plot(time_min, angVel(:,3) * rad2deg, 'b--');
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
