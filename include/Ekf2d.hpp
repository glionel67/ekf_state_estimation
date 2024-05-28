///
/// \file Ekf2d.hpp
/// \brief 2D state estimation using an Extended Kalman filter.
/// \author Lionel GENEVE
/// \date 17/11/2018
/// \version 1.0
///

#pragma once

// ---------- Headers ----------
// STD/C++ headers
#include <iostream>
#include <vector>

// EIGEN headers
#include <Eigen/Core>
#include <Eigen/Dense>

// Project headers
#include "Ekf.hpp"
#include "State2d.hpp"

namespace ekf {

class Ekf2d : public Ekf<State2d> {
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  Ekf2d() = delete;  ///< Default constructor

  Ekf2d(const FilterParameters& _params);  ///< Constructor

  Ekf2d(const Ekf2d& _ekf) = default;  ///< Default copy constructor

  virtual ~Ekf2d() = default;  ///< Default destructor

  void init(const VectorX& _x0, const MatrixX& _P0, double _tStamp);

  bool checkSanity(void);  ///< Check filter sanity

 protected:
  void predict(const State2d* _currState, State2d* _nextState);

  bool correctAccelerometer(State2d* _state, const Measurement* _meas);
  bool correctGyroscope(State2d* _state, const Measurement* _meas);
  bool correctAccGyr(State2d* _state, const Measurement* _meas);
  bool correctLinVel(State2d* _state, const Measurement* _meas);
  bool correctAngVel(State2d* _state, const Measurement* _meas);
  bool correctLinAngVel(State2d* _state, const Measurement* _meas);
  bool correctLinAcc(State2d* _state, const Measurement* _meas);
  bool correctPosition(State2d* _state, const Measurement* _meas);
  bool correctAttitude(State2d* _state, const Measurement* _meas);
  //bool correctGpsPosition(State2d* _state, const Measurement* _meas);
  //bool correctGpsSpeed(State2d* _state, const Measurement* _meas);

  void addErrorState(VectorX& _x, VectorX& _dx);
  void errorReset(MatrixX& _P, VectorX& _dx);
};  // class Ekf2d

}  // namespace ekf
