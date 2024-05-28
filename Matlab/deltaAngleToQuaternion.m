function q = deltaAngleToQuaternion(deltaAngle)
%deltaAngleToQuaternion: convert an angle variation into a quaternion form
% Args:
% - deltaAngle (vector [3x1]): axis/angle rotation
% Return:
% - q (vector [4x1]): quaternion form

rotAngle = 0.5 * norm(deltaAngle); % Rotation angle
rotAxis = deltaAngle; % Rotation axis

if (rotAngle >= 1e-5)
    rotAxis = deltaAngle / norm(deltaAngle);
end

ca = cos(rotAngle);
sa = sin(rotAngle);

q = [ca; rotAxis(1) * sa; rotAxis(2) * sa; rotAxis(3) * sa];

end
