/**
 * \file GpsConversion.cpp
 * \author Lionel GENEVE
 * \version 1.0
 * \date 17/11/2019
 * \brief GPS conversion functions
 * */

#include "Gps.hpp"
#include "GpsConversion.hpp"


void llaToEcef(const double _lat, const double _lon, const double _alt,
                double* _ecefX, double* _ecefY, double* _ecefZ)
{
    double sinLat = sin(_lat);
    double sinLat2 = sinLat * sinLat;
    double N = SEMI_MAJOR_AXIS / sqrt(1. - FIRST_ECCENTRICITY_SQUARED * sinLat2);
    double tmp = (N + _alt) * cos(_lat);
    
    *_ecefX = tmp * cos(_lon);
    *_ecefY = tmp * sin(_lon);
    *_ecefZ = (N * (1. - FIRST_ECCENTRICITY_SQUARED) + _alt) * sinLat;
}

void ecefToEnu(const double _ecefX, const double _ecefY, const double _ecefZ,
                const double _latRef, const double _lonRef, const double _altRef,
                double* _east, double* _north, double* _up)
{
    // Convert reference LLA to reference ECEF
    double ecefXRef = 0., ecefYRef = 0., ecefZRef = 0.;
    llaToEcef(_latRef, _lonRef, _altRef, &ecefXRef, &ecefYRef, &ecefZRef);

    double dx = _ecefX - ecefXRef;
    double dy = _ecefY - ecefYRef;
    double dz = _ecefZ - ecefZRef;

    double sinLat = sin(_latRef);
    double cosLat = cos(_latRef);

    double sinLon = sin(_lonRef);
    double cosLon = cos(_lonRef);

    *_east = -sinLon * dx + cosLon * dy;
    *_north = -sinLat * cosLon * dx - sinLat * sinLon * dy + cosLat * dz;
    *_up = cosLat * cosLon * dx + cosLat * sinLon * dy + sinLat * dz;
}
