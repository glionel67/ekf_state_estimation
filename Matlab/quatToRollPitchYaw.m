function [roll, pitch, yaw] = quatToRollPitchYaw(q)
%quatToRollPitchYaw: convert the quaternion <q> to roll/pitch/yaw angles in radians
% Args:
% - q (vector [4x1]): quaternion
% Return:
% - roll (double): roll angle [rad]
% - pitch (double): pitch angle [rad]
% - yaw (double): yaw angle [rad]

assert(length(q) == 4, 'quatToRollPitchYaw: q must be a [4x1] vector!');

% Extract quaternion components
qw = q(1);
qx = q(2);
qy = q(3);
qz = q(4);

% Compute roll/pitch/yaw angles [rad]
roll = atan2(2.0 * (qw * qx + qy * qz), 1.0 - 2.0 * (qx^2 + qy^2));
pitch = -asin(2.0 * (qw * qy - qz * qx));
yaw = atan2(2.0 * (qw * qz + qx * qy), 1.0 - 2.0 * (qy^2 + qz^2));

end
