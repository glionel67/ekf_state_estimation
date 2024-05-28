///
/// \file Ekf3d.hpp
/// \brief 3D state estimation using an Extended Kalman filter.
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
#include "State3d.hpp"

namespace ekf {

class Ekf3d : public Ekf<State3d> {
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  Ekf3d() = delete;  ///< Default constructor

  Ekf3d(const FilterParameters& _params);  ///< Constructor

  Ekf3d(const Ekf3d& _ekf) = default;  ///< Default copy constructor

  virtual ~Ekf3d() = default;  ///< Default destructor

  void init(const VectorX& _x0, const MatrixX& _P0, double _tStamp);

  bool checkSanity(void);  ///< Check filter sanity

 protected:
  void predict(const State3d* _currState, State3d* _nextState);

  bool correctAccelerometer(State3d* _state, const Measurement* _meas);
  bool correctGyroscope(State3d* _state, const Measurement* _meas);
  bool correctAccGyr(State3d* _state, const Measurement* _meas);
  bool correctLinVel(State3d* _state, const Measurement* _meas);
  bool correctAngVel(State3d* _state, const Measurement* _meas);
  bool correctLinAngVel(State3d* _state, const Measurement* _meas);
  bool correctLinAcc(State3d* _state, const Measurement* _meas);
  bool correctPosition(State3d* _state, const Measurement* _meas);
  bool correctAttitude(State3d* _state, const Measurement* _meas);
  //bool correctGpsPosition(State3d* _state, const Measurement* _meas);
  //bool correctGpsSpeed(State3d* _state, const Measurement* _meas);

  void addErrorState(VectorX& _x, VectorX& _dx);
  void errorReset(MatrixX& _P, VectorX& _dx);

  ///
  /// \fn deltaAngleToQuaternion
  /// \brief Convert a delta angle vector to a quaternion
  /// \param _deltaAngle: delta angle vector to convert (angles are in radians)
  /// \return Quaternion: the corresponding quaternion representation
  ///
  Quaternion deltaAngleToQuaternion(const Vector3& _deltaAngle);  // Convert a delta angle vector to a quaternion

  /**
   * @brief Convert a quaterion representation to Euler angles ZYX roll/pitch/yaw
   * @param _q quaterion to convert
   * @param _roll roll angle (rotation around x) [rad]
   * @param _pitch pitch angle (rotation around y) [rad]
   * @param _yaw yaw angle (rotation around z) [rad]
   */
  void quaternionToRollPitchYaw(const Quaternion& _q, double* _roll, double* _pitch, double* _yaw);
};  // class Ekf3d

}  // namespace ekf
