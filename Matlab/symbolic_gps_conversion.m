syms latRef lonRef altRef % lat=phi, lon=lambda
syms x y z

sinLat = sin(latRef);
cosLat = cos(latRef);

sinLon = sin(lonRef);
cosLon = cos(lonRef);

R = [-sinLat*cosLon, -sinLat*sinLon, cosLat;
    -sinLon, cosLon, 0;
    -cosLat*cosLon, -cosLat*sinLon, -sinLat];

ned = R * [x; y; z]