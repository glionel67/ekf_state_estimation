function [ecef] = lla2ecef(lla)
%lla2ecef: convert GPS coordinates to ECEF coordinates
% Args:
% - lla (matrix [3xn]): lat/lon/alt coordinates [rad/rad/m]
% Return:
% - ecef (matrix [3xn]): ECEF coordinates [m/m/m]

% Parameters
a = 6378137.0; % major radius of the earth reference ellipsoid (m)
f = 1.0 / 298.257223563; % flattening of the earth reference ellipsoid
e2 = f * (2 - f);

lat = lla(:,1);
lon = lla(:,2);
alt = lla(:,3);

v = a ./ sqrt(1.0 - e2 * sin(lat).^2);

tmp = (v + alt) .* cos(lat);
ecefX = tmp .* cos(lon);
ecefY = tmp .* sin(lon);
ecefZ = v * (1.0 - e2) .* sin(lat);

ecef = [ecefX, ecefY, ecefZ];
end