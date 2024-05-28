% EKF estimation filter symbolic computations
clear all;
close all;
clc;

nStates = 3 + 4 + 3 + 3 + 3 + 3 + 3
nErrorStates = 3 + 3 + 3 + 3 + 3 + 3 + 3

syms qw qx qy qz % Quaternion
syms gz % Gravity

R = [qw^2+qx^2-qy^2-qz^2, 2.0*(qx*qy-qw*qz), 2.0*(qx*qz+qw*qy);
    2.0*(qx*qy+qw*qz), qw^2-qx^2+qy^2-qz^2, 2.0*(qy*qz-qw*qx);
    2.0*(qx*qz-qw*qy), 2.0*(qy*qz+qw*qx), qw^2-qx^2-qy^2+qz^2]
g = [0.0; 0.0; gz]

% Measurement Jacobian matrix
Hdx = zeros(nStates, nErrorStates);
Hdx(1:3, 1:3) = eye(3);

Qdt = 0.5 * [-qx, -qy, -qz;
    qw, -qz, qy;
    qz, qw, -qx;
    -qy, qx, qw];
Hdx = sym(Hdx);
Hdx(4:7, 4:6) = Qdt;
Hdx(8:10, 7:9) = eye(3);
Hdx(11:13, 10:12) = eye(3);
Hdx(14:16, 13:15) = eye(3);
Hdx(17:19, 16:18) = eye(3);
Hdx(20:22, 19:21) = eye(3);
Hdx

%% Lin vel measurement jacobian matrix
Hx = zeros(3, nStates);
Hx(1:3, 8:10) = eye(3);
Hx = sym(Hx)
H_lin_vel = Hx * Hdx

%% Lin vel measurement jacobian matrix
Hx = zeros(3, nStates);
Hx(1:3, 11:13) = eye(3);
Hx = sym(Hx)
H_ang_vel = Hx * Hdx

%% Lin + ang vel measurement jacobian matrix
Hx = zeros(3, nStates);
Hx(1:3, 8:10) = eye(3);
Hx(4:6, 11:13) = eye(3);
Hx = sym(Hx)
H_lin_ang_vel = Hx * Hdx

%% Accelero measurement jacobian matrix
Hx = zeros(3, nStates);
Hx(1:3, 14:16) = eye(3);
Hx(1:3, 17:19) = eye(3);
Hx = sym(Hx);
Rtg = transpose(R) * g

Hx(1:3, 4:7) = [...
    diff(Rtg(1), qw), diff(Rtg(1), qx), diff(Rtg(1), qy), diff(Rtg(1), qz);
    diff(Rtg(2), qw), diff(Rtg(2), qx), diff(Rtg(2), qy), diff(Rtg(2), qz);
    diff(Rtg(3), qw), diff(Rtg(3), qx), diff(Rtg(3), qy), diff(Rtg(3), qz)]
H_accelero = Hx * Hdx
H_accelero_att = H_accelero(1:3, 4:6)

Rtgx = [0.0, -Rtg(3), Rtg(2);
    Rtg(3), 0.0, -Rtg(1);
    -Rtg(2), Rtg(1), 0.0]
% Check equality
simplify(H_accelero_att(1,1) - Rtgx(1,1))
simplify(H_accelero_att(1,2) - Rtgx(1,2))
simplify(H_accelero_att(1,3) - Rtgx(1,3))
simplify(H_accelero_att(2,1) - Rtgx(2,1))
simplify(H_accelero_att(2,3) - Rtgx(2,3))

%% Attitude/Yaw measurement
num = 2.0 * (qw * qx + qy * qz);
den = 1.0 - 2.0 * (qx^2 + qy^2);
roll_pred = atan(num /den)
pitch_pred = asin(2.0 * (qw * qy - qz * qx))
num = 2.0 * (qw * qz + qx * qy);
den = 1.0 - 2.0 * (qy^2 + qz^2);
yaw_pred = atan(num / den)
Hx = zeros(3, nStates);
Hx = sym(Hx);
Hx(1, 4:7) = [diff(roll_pred, qw), diff(roll_pred, qx), diff(roll_pred, qy), diff(roll_pred, qz)]
Hx(2, 4:7) = [diff(pitch_pred, qw), diff(pitch_pred, qx), diff(pitch_pred, qy), diff(pitch_pred, qz)]
Hx(3, 4:7) = [diff(yaw_pred, qw), diff(yaw_pred, qx), diff(yaw_pred, qy), diff(yaw_pred, qz)]
H_rpy = simplify(Hx * Hdx)

Hx = zeros(1, nStates);
Hx = sym(Hx);
Hx(1, 4:7) = [diff(yaw_pred, qw), diff(yaw_pred, qx), diff(yaw_pred, qy), diff(yaw_pred, qz)]
H_yaw = simplify(Hx * Hdx)


% Covariance matrix error reset
P = sym('P', [nErrorStates, nErrorStates]);

G = eye(nErrorStates, nErrorStates);
G = sym(G)
Ga = sym('Ga', [3, 3])
G(4:6,4:6) = Ga

GPGt = G * P * transpose(G)

P == GPGt


%% Prediction step
nProcNoises = 12; % 3: ang vel, 3: lin acc, 3: acc bias, 3: gyr bias
W = zeros(nErrorStates, nProcNoises);
W(10:12, 1:3) = eye(3);
W(13:15, 4:6) = eye(3);
W(16:18, 7:9) = eye(3);
W(19:21, 10:12) = eye(3)

Qv = sym('Qv', [nProcNoises, 1]);
Q = diag(Qv)

WQWt = W * Q * transpose(W)