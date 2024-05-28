function R = rollPitchYawToRotMat(roll, pitch, yaw)
%rollPitchYawToRotMat: convert the roll/pitch/yaw Euler angles in [rad] to
%a rotation matrix (Direct Cosinus Matrix or DCM)

Rz = [cos(yaw), -sin(yaw), 0.0;
    sin(yaw), cos(yaw), 0.0;
    0.0, 0.0, 1.0];

Ry = [cos(pitch), 0.0, sin(pitch);
    0.0, 1.0, 0.0;
    -sin(pitch), 0.0, cos(pitch)];
% Ry = [cos(pitch), 0.0, -sin(pitch);
%     0.0, 1.0, 0.0;
%     sin(pitch), 0.0, cos(pitch)]; ???????????????

Rx = [1.0, 0.0, 0.0;
    0.0, cos(roll), -sin(roll);
    0.0, sin(roll), cos(roll)];
% Rx = [1.0, 0.0, 0.0;
%     0.0, cos(roll), sin(roll);
%     0.0, -sin(roll), cos(roll)]; ???????????????

R = Rz * Ry * Rx;

end
