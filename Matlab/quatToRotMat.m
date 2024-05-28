function R = quatToRotMat(q)
%quatToRotMat: convert a quaternion q to a rotation matrix R
% Args:
% - q (vector [4x1]): quaternion
% Return:
% - R (matrix [3x3]): rotation matrix \in SO(3)

assert(length(q) == 4, 'quatToRotMat: quaternion <q> must be a [4x1] vector!');

% First row
r11 = q(1)^2 + q(2)^2 - q(3)^2 - q(4)^2;
r12 = 2*(q(2)*q(3) - q(1)*q(4));
r13 = 2*(q(2)*q(4) + q(1)*q(3));
% Second row
r21 = 2*(q(2)*q(3) + q(1)*q(4));
r22 = q(1)^2 - q(2)^2 + q(3)^2 - q(4)^2;
r23 = 2*(q(3)*q(4) - q(1)*q(2));
% Third row
r31 = 2*(q(2)*q(4) - q(1)*q(3));
r32 = 2*(q(3)*q(4) + q(1)*q(2));
r33 = q(1)^2 - q(2)^2 - q(3)^2 + q(4)^2;
% Rotation matrix
R = [r11, r12, r13;
    r21, r22, r23;
    r31, r32, r33];

end
