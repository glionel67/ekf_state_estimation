///
/// \file FilterParameters.hpp
/// \brief Estimation filter (EKF) parameters
/// \author Lionel GENEVE
/// \date 17/11/2018
/// \version 1.0
///

#pragma once

// ---------- Headers ----------
// STD/C++ headers
//#include <iostream>
#include <string>

// EIGEN headers
#include <Eigen/Dense>

// Project headers
#include "FilterUtils.hpp"

///
/// \namespace ekf
/// \brief ekf namespace
///
namespace ekf {

///
/// \class FilterParameters
/// \brief Parameters of the filter
///
class FilterParameters {
public:
  FilterParameters() = default; ///< Default constructor
  ~FilterParameters() = default; ///< Default destructor
  FilterParameters(const FilterParameters& _filterParams) = default; ///< Default copy constructor
  FilterParameters& operator=(const FilterParameters& _filterParams) = default;

  bool loadParameters(const std::string& _filename); // TODO: load parameters from YAML file

  bool checkParameters(void); // TODO: check parameters validity
  
  // Parameters
  double positionProcessNoiseStd;
  double attitudeProcessNoiseStd;
  double linVelProcessNoiseStd;
  double linAccProcessNoiseStd;
  double angVelProcessNoiseStd;
  double accBiasNoiseStd;  ///< Accelerometer bias noise standard deviation
  double gyrBiasNoiseStd;  ///< Gyroscope bias noise standard deviation

  size_t nIters;  ///< Number of iterations for the correction step (iterated EKF)
  double temporalWindow;  ///< Temporal window of the state buffer [s]
  double minCovMatrix, maxCovMatrix;  ///< Min/max covariance matrix bounds
  Vector3 minPosition, maxPosition;  ///< 3D position bounds [m]
  Vector3 minAttitude, maxAttitude;  ///< 3D attitude bounds [rad]
  Vector3 minLinVel, maxLinVel;  ///< 3D linear velocity bounds [m/s]
  Vector3 minAngVel, maxAngVel;  ///< 3D angular velocity bounds [rad/s]
  Vector3 minLinAcc, maxLinAcc;  ///< 3D linear acceleration bounds [m/s2]
  Vector3 minAccBias, maxAccBias; ///< Accelerometer bias bounds [m/s2]
  Vector3 minGyrBias, maxGyrBias; ///< Gyroscope bias bounds [rad/s]

  Vector3 maxPositionStd;  ///< Maximal 3D position estimate standard deviation
  Vector3 maxAttitudeStd;  ///< Maximal 3D attitude estimate standard deviation
  Vector3 maxLinVelStd;  ///< Maximal 3D linear velocity estimate standard deviation
  Vector3 maxAngVelStd;  ///< Maximal 3D angular velocity estimate standard deviation
  Vector3 maxLinAccStd;  ///< Maximal 3D linear acceleration estimate standard deviation
  Vector3 maxAccBiasStd;  ///< Maximal accelero bias estimate standard deviation
  Vector3 maxGyrBiasStd;  ///< Maximal gyro bias estimate standard deviation

  Vector3 gravityVector; ///< Earth gravity vector (normally [0; 0; -9.81])
};  // class FilterParameters

} // namespace ekf
