function angle = wrapAngle(angle)
%wrapAngle: keep the given angle in (-pi, pi]
% Args:
% - angle (double or vector): angle to wrap [rad]
% Return:
% - angle (double or vector): wrapped angle in (-pi, pi] [rad]

% Remove multiples of 2*pi
% a = mod(abs(angle), 2.0 * pi);
% 
% mask = angle < 0.0;
% angle(mask) = - angle(mask);
% if angle < 0.0
%     angle = - angle;
% end

% Keep angle in (-pi, pi]
mask = angle > pi;
angle(mask) = angle(mask) - 2.0 * pi;
mask = angle <= -pi;
angle(mask) = angle(mask) + 2.0 * pi;

% if (angle > pi)
%     angle = angle - 2.0 * pi;
% elseif (angle <= -pi)
%     angle = angle + 2.0 * pi;
% %else
%     %nop;
% end

end
