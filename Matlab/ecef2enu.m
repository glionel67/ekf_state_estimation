function [enu] = ecef2enu(ecef, llaRef)
%ecef2enu: convert ECEF (Earth Center Earth Fixed) coordinates to ENU (East North Up) coordinates
% Args:
% - ecef (matrix [3xn]): ECEF coordinates [m/m/m]
% - llaRef (vector [3x1]): lat/lon/alt reference for local position [rad/rad/m]
% Return:
% - enu (matrix [3xn]): ENU coordinates [m/m/m]

sLat = sin(llaRef(1));
cLat = cos(llaRef(1));
sLon = sin(llaRef(2));
cLon = cos(llaRef(2));

R = [-sLon, cLon, 0;
    -sLat * cLon, -sLat * sLon, cLat;
    cLat * cLon, cLat * sLon, sLat];

%disp(['R=', num2str(R)]);

ecefRef = lla2ecef(llaRef);

enu = R * (ecef' - ecefRef');

enu = enu';
end

