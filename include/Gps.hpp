/**
 * \file Gps.hpp
 * \author Lionel GENEVE
 * \version 1.0
 * \date 17/11/2019
 * \brief GPS useful constants (based on WGS-84 model)
 * */
#pragma once

#include <math.h>

// Angle conversions
#ifndef DEG_TO_RAD
#define DEG_TO_RAD (M_PI / 180.)
#endif
#ifndef RAD_TO_DEG
#define RAD_TO_DEG (180. / M_PI)
#endif

#define EARTH_ROTATION_RATE (7.2921151467e-05) // Earth rotation rate [rad/s]

#define EARTH_GRAVITY (9.812) // Earth gravity norm [m/s^2]

#define LIGHT_SPEED (299792458.0) // Speed of light [m/s]

#define SEMI_MAJOR_AXIS (6378137.) // WGS84 equatorial earth radius/major radius of the earth reference ellipsoid [m] (WGS84 datum)
#define SEMI_MINOR_AXIS (6356752.314245) // [m]
#define FIRST_ECCENTRICITY_SQUARED (.00669437999014132) // First eccentricity squared: e^2 = f*(2-f)
#define SECOND_ECCENTRICITY_SQUARED (.00673949674228)
#define FLATTENING (1./298.257223563) // flattening of the earth reference ellipsoid  (WGS84 datum)
