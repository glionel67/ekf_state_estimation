///
/// \file FilterUtils.cpp
/// \brief Filter utility defines and functions
/// \author Lionel GENEVE
/// \date 17/11/2018
/// \version 1.0
///

// ---------- Headers ----------
// STD/C++ headers
#include <iostream>

// EIGEN headers

// Project headers
#include "FilterUtils.hpp"

///
/// \namespace ekf
/// \brief ekf namespace
///
namespace ekf {

inline void constrainValue(double& _val, const double& _min, const double& _max) {
    if (std::isnan(_val) || std::isinf(_val)) {
        std::cerr << "constrainValue: isnan or isinf detected!\n";
        _val  = 0.;
    }
    else if (_val > _max) {
        _val = _max;
    }
    else if (_val < _min) {
        _val = _min;
    }
    //else { }
}  // constrainValue

double piToPi(double _angle) {
  // Remove multiples of 2*pi
  double angle = std::fmod(std::fabs(_angle), 2.0 * M_PI);
  if (_angle < 0.0) angle = -angle;

  // Keep angle in (-pi, pi]
  if (angle > M_PI) {
    angle -= 2.0 * M_PI;
  } else if (angle <= - M_PI) {
    angle += 2.0 * M_PI;
  }

  return angle;
}  // piToPi

} // namespace ekf