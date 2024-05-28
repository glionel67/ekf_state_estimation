% Create EKF dataset

clear all;
close all;
clc;

%% Parameters
steerRawPosToMm = 0.008;
deg2rad = pi / 180.0;
rad2deg = 180.0 / pi;
ms2kmh = 3.6;
kmh2ms = 1.0 / 3.6;

%% Select data
path = 'D:\HE441-Data\2022-03-24';
%filename = '20220324_100451.mat';
%filename = '20220324_134720.mat';
filename = '20220324_134720_cleaned.mat';

%% Load data
load(fullfile(path, filename));

if ~exist('Time', 'var')
    error('ERROR: no Time signal');
end
if ~exist('Time_bis', 'var') % For CAN 5 with CSM only
    warning('Creating fake --Time_bis-- variable');
    Time_bis = Time;
end

nData = length(Time);
timeMin = Time / 60; % Convert time in [s] to [min]

%% Select data range
rng = 1:nData; % Select all data
% Test on 01/03/2022
%rng = (timeMin > 40) & (timeMin < 50);
%rng = (timeMin > 80) & (timeMin < 100);
% Test on 24/03/2022
%rng = (timeMin > 26) & (timeMin < 40);
%rng = (timeMin > 104) & (timeMin < 110);
%rng = (timeMin > 118) & (timeMin < 128);

%% Plot

figure;
subplot(3,1,1);
plot(Time_bis(rng)/60, VEHICLE_ANGULAR_VELOCITY_VehicleRollRate(rng), 'r-'); hold on;
plot(Time_bis(rng)/60, RateVehicle_AngRateX(rng)*deg2rad, 'b--');
plot(Time_bis(rng)/60, RateLevel_AngRateForward(rng)*deg2rad, 'g:');
legend('PX4', 'RT2000');
ylabel('Ang vel X [rad/s]');
title('IMUs gyroscope comparison');
subplot(3,1,2);
plot(Time_bis(rng)/60, VEHICLE_ANGULAR_VELOCITY_VehiclePitchRate(rng), 'r-'); hold on;
plot(Time_bis(rng)/60, RateVehicle_AngRateY(rng)*deg2rad, 'b--');
plot(Time_bis(rng)/60, RateLevel_AngRateLateral(rng)*deg2rad, 'g:');
legend('PX4', 'RT2000');
ylabel('Ang vel Y [rad/s]');
subplot(3,1,3);
plot(Time_bis(rng)/60, VEHICLE_ANGULAR_VELOCITY_VehicleYawRate(rng), 'r-'); hold on;
plot(Time_bis(rng)/60, RateVehicle_AngRateZ(rng)*deg2rad, 'b--');
plot(Time_bis(rng)/60, RateLevel_AngRateDown(rng)*deg2rad, 'g:');
legend('PX4', 'RT2000');
ylabel('Ang vel Z [rad/s]');
xlabel('Time [s]');

figure;
subplot(3,1,1);
plot(Time_bis(rng)/60.0, VEHICLE_LINEAR_ACCEL_VehicleLinearAccelX(rng), 'r-'); hold on;
plot(Time_bis(rng)/60.0, AccelVehicle_AccelX(rng), 'b--');
plot(Time_bis(rng)/60.0, AccelLevel_AccelForward(rng), 'g:');
legend('PX4', 'RT2000');
ylabel('Acc X [m/s2]');
title('IMUs accelerometer comparison');
subplot(3,1,2);
plot(Time_bis(rng)/60.0, VEHICLE_LINEAR_ACCEL_VehicleLinearAccelY(rng), 'r-'); hold on;
plot(Time_bis(rng)/60.0, AccelVehicle_AccelY(rng), 'b--');
plot(Time_bis(rng)/60.0, AccelLevel_AccelLateral(rng), 'g:');
legend('PX4', 'RT2000');
ylabel('Acc Y [m/s2]');
subplot(3,1,3);
plot(Time_bis(rng)/60.0, VEHICLE_LINEAR_ACCEL_VehicleLinearAccelZ(rng), 'r-.');  hold on;
plot(Time_bis(rng)/60.0, AccelVehicle_AccelZ(rng), 'b--');
plot(Time_bis(rng)/60.0, AccelLevel_AccelDown(rng), 'g:');
legend('PX4', 'RT2000');
ylabel('Acc Z [m/s2]');
xlabel('Time [min]');

figure;
plot(Time_bis(rng)/60.0, Velocity_Speed2D(rng)*3.6, 'b-'); hold on;
plot(Time(rng)/60.0, VEHICLE_STATE_1_VehicleSpeed(rng), 'r--');
legend('Speed2D', 'Wheel based');
ylabel('Velocity [km/h]');
xlabel('Time [min]');

figure; hold on;
plot(Time_bis(rng)/60.0, MOTOR_SPEEDS_MotorSpeedFL(rng), 'r-');
plot(Time_bis(rng)/60.0, MOTOR_SPEEDS_MotorSpeedFR(rng), 'b--');
plot(Time_bis(rng)/60.0, MOTOR_SPEEDS_MotorSpeedRL(rng), 'g:');
plot(Time_bis(rng)/60.0, MOTOR_SPEEDS_MotorSpeedRR(rng), 'c-.');
legend('FL', 'FR', 'RL', 'RR');
ylabel('Wheel speed [rpm]');
xlabel('Time [min]');
title('Vehicle motor/wheel speeds');

figure; hold on;
plot(Time_bis(rng)/60.0, VEHICLE_DIRECTION_ANGLES_VehicleDirectionAngleFront(rng), 'r-');
plot(Time_bis(rng)/60.0, VEHICLE_DIRECTION_ANGLES_VehicleDirectionAngleRear(rng), 'b--');
legend('Front', 'Rear');
ylabel('Steering angle [rad]');
xlabel('Time [min]');
title('Vehicle steering angles');

figure;
plot(GPS_LAT_LON_Longitude(rng), GPS_LAT_LON_Latitude(rng), 'b-'); hold on;
plot(LatitudeLongitude_PosLon(rng), LatitudeLongitude_PosLat(rng), 'g--');
legend('PX4', 'RT2000');
xlabel('Longitude [deg]');
ylabel('Latitude [deg]');
title('GPS trajectory');

figure;
plot(Time_bis(rng)/60.0, GPS_ALT_FIX_SAT_Altitude(rng), 'b-'); hold on;
plot(Time_bis(rng)/60.0, Altitude_PosAlt(rng), 'r--');
legend('PX4', 'RT2000');
xlabel('Time [min]');
ylabel('Altitude [m]');
title('GPS altitude');

figure;
plot(Time_bis(rng)/60.0, GPS_ALT_FIX_SAT_FixType(rng), 'b+'); hold on;
plot(Time_bis(rng)/60.0, GpsStatus_GpsPosMode(rng), 'ro');
legend('PX4', 'RT2000');
xlabel('Time [min]');
ylabel('Fix type');
title('GPS fix type/mode');

%% Create dataset

% Get measurements
% RT2000 IMU measurements
zAcc1 = [AccelVehicle_AccelX(rng), ...
    AccelVehicle_AccelY(rng), ...
    AccelVehicle_AccelZ(rng)]; % RT2000
zGyr1 = [RateVehicle_AngRateX(rng)*deg2rad, ...
    RateVehicle_AngRateY(rng)*deg2rad, ...
    RateVehicle_AngRateZ(rng)*deg2rad]; % RT2000
% PX4 IMU measurements
zAcc2 = [VEHICLE_LINEAR_ACCEL_VehicleLinearAccelX(rng), ...
    VEHICLE_LINEAR_ACCEL_VehicleLinearAccelY(rng), ...
    VEHICLE_LINEAR_ACCEL_VehicleLinearAccelZ(rng)]; % PX4
zGyr2 = [VEHICLE_ANGULAR_VELOCITY_VehicleRollRate(rng), ...
    VEHICLE_ANGULAR_VELOCITY_VehiclePitchRate(rng), ...
    VEHICLE_ANGULAR_VELOCITY_VehicleYawRate(rng)]; % PX4
% Wheel motor speeds [rpm]
zWheels = [MOTOR_SPEEDS_MotorSpeedFL(rng), ...
    MOTOR_SPEEDS_MotorSpeedFR(rng), ...
    MOTOR_SPEEDS_MotorSpeedRL(rng), ...
    MOTOR_SPEEDS_MotorSpeedRR(rng)]; % [rpm]
% Steering angles [rad]
zSteers = [VEHICLE_DIRECTION_ANGLES_VehicleDirectionAngleFront(rng), ...
    VEHICLE_DIRECTION_ANGLES_VehicleDirectionAngleRear(rng)]; % [rad]
% GPS measurements
zGps = [GPS_LAT_LON_Latitude(rng), ...
    GPS_LAT_LON_Longitude(rng), ...
    GPS_ALT_FIX_SAT_Altitude(rng), ...
    GPS_ALT_FIX_SAT_FixType(rng), ...
    GPS_ALT_FIX_SAT_NumberOfSatellites(rng)];
% Vehicle measurements
zVeh = [VEHICLE_STATE_1_TransmissionMode(rng), ... % 0: sna, 1: park, 2: reverse, 3: neutral, 4: forward
    VEHICLE_STATE_1_VehicleSpeed(rng)];
% Create dataset
dataset = [Time(rng), zAcc1, zGyr1, zAcc2, zGyr2, zWheels, zSteers, zGps, zVeh];
% Save dataset
f = 'ekfDataset4.csv';
writematrix(dataset, f, 'FileType', 'text', 'Delimiter', ';');


%% Plot results of state estimation EKF filter
%f = 'D:\Workspace\StateEstimationFilter\Data\ekfResults.csv';
%f = 'D:\Workspace\StateEstimationFilter\Data\ekfResults0.csv';
f = 'D:\Workspace\StateEstimationFilter\Data\ekfResults2b.csv';
res = readmatrix(f);
p_ekf = res(:, 1:3);
q_ekf = res(:, 4:7);
v_ekf = res(:, 8:10);
w_ekf = res(:, 11:13);
a_ekf = res(:, 14:16);
ba_ekf = res(:, 17:19);
bg_ekf = res(:, 20:22);
rpy_ekf = zeros(size(p_ekf));
for i=1:size(rpy_ekf, 1)
    rpy_ekf(i,:) = quatToRollPitchYaw(q_ekf(i,:));
end

p_std = sqrt(res(:, 23:25));
q_std = sqrt(res(:, 26:28));
v_std = sqrt(res(:, 29:31));
w_std = sqrt(res(:, 32:34));
a_std = sqrt(res(:, 35:37));
ba_std = sqrt(res(:, 38:40));
bg_std = sqrt(res(:, 41:43));

lla = [LatitudeLongitude_PosLat(rng) * deg2rad,...
    LatitudeLongitude_PosLon(rng) * deg2rad,...
    Altitude_PosAlt(rng) * deg2rad];
ecef = lla2ecef(lla);
enu = ecef2enu(ecef, lla(1,:));

figure;
plot(enu(:,1), enu(:,2), 'r-'); hold on;
plot(p_ekf(:,1), p_ekf(:,2), 'b--'); hold on;
plot(p_ekf(1,1), p_ekf(1,2), 'ks');
plot(p_ekf(end,1), p_ekf(end,2), 'gd');
legend('RT2000', 'EKF');
xlabel('X [m]');
ylabel('Y [m]');
title('Position');

figure;
subplot(3,1,1);
plot(Time_bis(rng)/60.0, HeadingPitchRoll_AngleRoll(rng), 'r-'); hold on;
plot(Time_bis(rng)/60.0, rpy_ekf(:,1) * rad2deg, 'b--');
legend('RT2000', 'EKF');
ylabel('Roll [deg]');
title('Attitude');
subplot(3,1,2);
plot(Time_bis(rng)/60.0, HeadingPitchRoll_AnglePitch(rng), 'r-'); hold on;
plot(Time_bis(rng)/60.0, rpy_ekf(:,2) * rad2deg, 'b--');
legend('RT2000', 'EKF');
ylabel('Pitch [deg]');
subplot(3,1,3);
plot(Time_bis(rng)/60.0, HeadingPitchRoll_AngleHeading(rng), 'r-.');  hold on;
plot(Time_bis(rng)/60.0, rpy_ekf(:,3) * rad2deg, 'b--');
legend('RT2000', 'EKF');
ylabel('Yaw [deg]');
xlabel('Time [min]');

figure;
subplot(3,1,1);
plot(Time_bis(rng)/60.0, VelocityLevel_VelForward(rng), 'r-'); hold on;
plot(Time_bis(rng)/60.0, v_ekf(:,1), 'b--');
legend('RT2000', 'EKF');
ylabel('Vel X [m/s]');
title('Linear velocity');
subplot(3,1,2);
plot(Time_bis(rng)/60.0, VelocityLevel_VelLateral(rng), 'r-'); hold on;
plot(Time_bis(rng)/60.0, v_ekf(:,2), 'b--');
legend('RT2000', 'EKF');
ylabel('Vel Y [m/s]');
subplot(3,1,3);
plot(Time_bis(rng)/60.0, Velocity_VelDown(rng), 'r-.');  hold on;
plot(Time_bis(rng)/60.0, v_ekf(:,3), 'b--');
legend('RT2000', 'EKF');
ylabel('Vel Z [m/s]');
xlabel('Time [min]');

figure;
subplot(3,1,1);
plot(Time_bis(rng)/60.0, RateLevel_AngRateForward(rng), 'r-'); hold on;
plot(Time_bis(rng)/60.0, w_ekf(:,1) * rad2deg, 'b--');
legend('RT2000', 'EKF');
ylabel('Vel X [deg/s]');
title('Angular velocity');
subplot(3,1,2);
plot(Time_bis(rng)/60.0, RateLevel_AngRateLateral(rng), 'r-'); hold on;
plot(Time_bis(rng)/60.0, w_ekf(:,2) * rad2deg, 'b--');
legend('RT2000', 'EKF');
ylabel('Vel Y [deg/s]');
subplot(3,1,3);
plot(Time_bis(rng)/60.0, RateLevel_AngRateDown(rng), 'r-.');  hold on;
plot(Time_bis(rng)/60.0, w_ekf(:,3) * rad2deg, 'b--');
legend('RT2000', 'EKF');
ylabel('Vel Z [deg/s]');
xlabel('Time [min]');

figure;
subplot(3,1,1);
plot(Time_bis(rng)/60.0, AccelLevel_AccelForward(rng), 'r-'); hold on;
plot(Time_bis(rng)/60.0, a_ekf(:,1), 'b--');
legend('RT2000', 'EKF');
ylabel('Acc X [m/s2]');
title('Linear acceleration');
subplot(3,1,2);
plot(Time_bis(rng)/60.0, AccelLevel_AccelLateral(rng), 'r-'); hold on;
plot(Time_bis(rng)/60.0, a_ekf(:,2), 'b--');
legend('RT2000', 'EKF');
ylabel('Acc Y [m/s2]');
subplot(3,1,3);
plot(Time_bis(rng)/60.0, AccelLevel_AccelDown(rng), 'r-.');  hold on;
plot(Time_bis(rng)/60.0, a_ekf(:,3), 'b--');
legend('RT2000', 'EKF');
ylabel('Acc Z [m/s2]');
xlabel('Time [min]');

figure;
subplot(3,1,1);
plot(Time_bis(rng)/60.0, ba_ekf(:,1), 'b--');
legend('EKF');
ylabel('X [m/s2]');
title('Accelero bias');
subplot(3,1,2);
plot(Time_bis(rng)/60.0, ba_ekf(:,2), 'b--');
legend('EKF');
ylabel('Y [m/s2]');
subplot(3,1,3);
plot(Time_bis(rng)/60.0, ba_ekf(:,3), 'b--');
legend('EKF');
ylabel('Z [m/s2]');
xlabel('Time [min]');

figure;
subplot(3,1,1);
plot(Time_bis(rng)/60.0, bg_ekf(:,1), 'b--');
legend('EKF');
ylabel('X [rad/s]');
title('Gyro bias');
subplot(3,1,2);
plot(Time_bis(rng)/60.0, bg_ekf(:,2), 'b--');
legend('EKF');
ylabel('Y [rad/s]');
subplot(3,1,3);
plot(Time_bis(rng)/60.0, bg_ekf(:,3), 'b--');
legend('EKF');
ylabel('Z [rad/s]');
xlabel('Time [min]');

figure; hold on;
plot(Time_bis(rng)/60.0, p_std(:,1), 'r-');
plot(Time_bis(rng)/60.0, p_std(:,2), 'b--');
plot(Time_bis(rng)/60.0, p_std(:,3), 'g:');
legend('x', 'y', 'z');
xlabel('Time [min]');
ylabel('Pos std [m]');
title('Position filter std');

figure; hold on;
plot(Time_bis(rng)/60.0, q_std(:,1), 'r-');
plot(Time_bis(rng)/60.0, q_std(:,2), 'b--');
plot(Time_bis(rng)/60.0, q_std(:,3), 'g:');
legend('x', 'y', 'z');
xlabel('Time [min]');
ylabel('Att std [rad]');
title('Attitude filter std');

figure; hold on;
plot(Time_bis(rng)/60.0, v_std(:,1), 'r-');
plot(Time_bis(rng)/60.0, v_std(:,2), 'b--');
plot(Time_bis(rng)/60.0, v_std(:,3), 'g:');
legend('x', 'y', 'z');
xlabel('Time [min]');
ylabel('Lin vel std [m/s]');
title('Linear velocity filter std');

figure; hold on;
plot(Time_bis(rng)/60.0, w_std(:,1), 'r-');
plot(Time_bis(rng)/60.0, w_std(:,2), 'b--');
plot(Time_bis(rng)/60.0, w_std(:,3), 'g:');
legend('x', 'y', 'z');
xlabel('Time [min]');
ylabel('Ang vel std [rad/s]');
title('Angular velocity filter std');

figure; hold on;
plot(Time_bis(rng)/60.0, a_std(:,1), 'r-');
plot(Time_bis(rng)/60.0, a_std(:,2), 'b--');
plot(Time_bis(rng)/60.0, a_std(:,3), 'g:');
legend('x', 'y', 'z');
xlabel('Time [min]');
ylabel('Lin acc std [m/s2]');
title('Linear acceleration filter std');

figure; hold on;
plot(Time_bis(rng)/60.0, ba_std(:,1), 'r-');
plot(Time_bis(rng)/60.0, ba_std(:,2), 'b--');
plot(Time_bis(rng)/60.0, ba_std(:,3), 'g:');
legend('x', 'y', 'z');
xlabel('Time [min]');
ylabel('Acc bias std [m/s2]');
title('Accelero bias filter std');

figure; hold on;
plot(Time_bis(rng)/60.0, bg_std(:,1), 'r-');
plot(Time_bis(rng)/60.0, bg_std(:,2), 'b--');
plot(Time_bis(rng)/60.0, bg_std(:,3), 'g:');
legend('x', 'y', 'z');
xlabel('Time [min]');
ylabel('Gyr bias std [rad/s]');
title('Gyro bias filter std');


%% Compute noise statistics
% mask=(Time_bis/60.0>47 & Time_bis/60.0<49.0);
% fprintf('----- RT2000 -----\n');
% % Mean <=> bias
% bax = mean(AccelVehicle_AccelX(mask));
% bay = mean(AccelVehicle_AccelY(mask));
% baz = mean(AccelVehicle_AccelZ(mask)) + 9.8;
% bgx = mean(RateVehicle_AngRateX(mask) * deg2rad);
% bgy = mean(RateVehicle_AngRateY(mask) * deg2rad);
% bgz = mean(RateVehicle_AngRateZ(mask) * deg2rad);
% fprintf('bax=%3.9f, bay=%3.9f, baz=%3.9f m/s2\n', bax, bay, baz);
% fprintf('bgx=%3.9f, bgy=%3.9f, bgz=%3.9f rad/s\n', bgx, bgy, bgz);
% % Std <=> std
% sax = std(AccelVehicle_AccelX(mask));
% say = std(AccelVehicle_AccelY(mask));
% saz = std(AccelVehicle_AccelZ(mask));
% sgx = std(RateVehicle_AngRateX(mask) * deg2rad);
% sgy = std(RateVehicle_AngRateY(mask) * deg2rad);
% sgz = std(RateVehicle_AngRateZ(mask) * deg2rad);
% fprintf('sax=%3.9f, say=%3.9f, saz=%3.9f m/s2\n', sax, say, saz);
% fprintf('sgx=%3.9f, sgy=%3.9f, sgz=%3.9f rad/s\n', sgx, sgy, sgz);
% fprintf('----- PX4 -----\n');
% % Mean <=> bias
% bax = mean(VEHICLE_LINEAR_ACCEL_VehicleLinearAccelX(mask));
% bay = mean(VEHICLE_LINEAR_ACCEL_VehicleLinearAccelY(mask));
% baz = mean(VEHICLE_LINEAR_ACCEL_VehicleLinearAccelZ(mask)) + 9.8;
% bgx = mean(VEHICLE_ANGULAR_VELOCITY_VehicleRollRate(mask) * deg2rad);
% bgy = mean(VEHICLE_ANGULAR_VELOCITY_VehiclePitchRate(mask) * deg2rad);
% bgz = mean(VEHICLE_ANGULAR_VELOCITY_VehicleYawRate(mask) * deg2rad);
% fprintf('bax=%3.9f, bay=%3.9f, baz=%3.9f m/s2\n', bax, bay, baz);
% fprintf('bgx=%3.9f, bgy=%3.9f, bgz=%3.9f rad/s\n', bgx, bgy, bgz);
% % Std <=> std
% sax = std(VEHICLE_LINEAR_ACCEL_VehicleLinearAccelX(mask));
% say = std(VEHICLE_LINEAR_ACCEL_VehicleLinearAccelY(mask));
% saz = std(VEHICLE_LINEAR_ACCEL_VehicleLinearAccelZ(mask));
% sgx = std(VEHICLE_ANGULAR_VELOCITY_VehicleRollRate(mask) * deg2rad);
% sgy = std(VEHICLE_ANGULAR_VELOCITY_VehiclePitchRate(mask) * deg2rad);
% sgz = std(VEHICLE_ANGULAR_VELOCITY_VehicleYawRate(mask) * deg2rad);
% fprintf('sax=%3.9f, say=%3.9f, saz=%3.9f m/s2\n', sax, say, saz);
% fprintf('sgx=%3.9f, sgy=%3.9f, sgz=%3.9f rad/s\n', sgx, sgy, sgz);

