% Script used to test the EKF implementation

clear all;
close all;
clc;

ts1 = 0.1; % Timestamp of the measurement [s]
z1 = [1; 2; 3]; % Measurement vector
R1 = diag([0.1; 0.1; 0.1]); % Measurement covariance matrix
id1 = 1; % Measurement unique ID
sensorName1 = 'DummySensor'; % Name of the corresponding sensor
mahaThresh1 = 8.0; % Mahalanobis threshold
m1 = Measurement(ts1, z1, R1, id1, sensorName1, mahaThresh1);
disp(m1.toString());

assert(isa(m1, 'Measurement'), 'Wrong class!');

m2 = Measurement(0.3, 2*z1, 0.1*R1, 3, 'UselessSensor', 3.99);
disp(m2.toString());

disp('m1 > m2 ?');
disp(m1 > m2);

disp('m1 < m2 ?');
disp(m1 < m2);

disp('m1 == m2 ?');
disp(m1 == m2);

ts1 = 0.0;
x1 = randn(10, 1);
P1 = randn(10, 10);
s1 = State(ts1, x1, P1);
disp(s1.toString());
assert(isa(s1, 'State'), 'Wrong class!');

s1.addMeasurement(m1);
disp(s1.toString());
s1.addMeasurement(m2);
disp(s1.toString());