/**
 * \file GpsConversion.hpp
 * \author Lionel GENEVE
 * \version 1.0
 * \date 17/11/2019
 * \brief GPS conversion functions
 * */
#pragma once

/**
 * \fn llaToEcef
 * \brief Convert lattitude/longitude/altitude to Earth-centered Earth-fixed coordinates
 * \param[in] _lat lattitude in radians
 * \param[in] _lon longitude in radians
 * \param[in] _alt altitude in meters
 * \param[out] _ecefX in meters
 * \param[out] _ecefY in meters
 * \param[out] _ecefZ in meters
 * \return none
 * */
void llaToEcef(const double _lat, const double _lon, const double _alt,
                double* _ecefX, double* _ecefY, double* _ecefZ);

/**
 * \fn ecefToEnu
 * \brief Convert Earth-centered Earth-fixed to East/North/Up coordinates
 * \param[in] _ecefX in meters [m]
 * \param[in] _ecefY in meters [m]
 * \param[in] _ecefZ in meters [m]
 * \param[in] _latRef reference lattitude in radians [rad]
 * \param[in] _lonRef reference longitude in radians [rad]
 * \param[in] _altRef reference altitude in meters [m]
 * \param[out] _east in meters [m]
 * \param[out] _north in meters [m]
 * \param[out] _up in meters [m]
 * \return none
 * */
void ecefToEnu(const double _ecefX, const double _ecefY, const double _ecefZ,
                const double _latRef, const double _lonRef, const double _altRef,
                double* _east, double* _north, double* _up);
