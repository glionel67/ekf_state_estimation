///
/// \file State3d.hpp
/// \brief Base class for an EKF state in 3D
/// \author Lionel GENEVE
/// \date 17/11/2018
/// \version 1.0
///
#pragma once

// ---------- Headers ----------
// STD/C++
#include <iostream>
#include <vector>
#include <string>
#include <memory>

// EIGEN
#include <Eigen/Dense>

// Project
#include "State.hpp"

// ---------- Defines ----------

///
/// \namespace ekf
/// \brief ekf namespace
///
namespace ekf {

class State3d : public State {
 public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    typedef std::unique_ptr<State3d> Ptr;
    typedef std::unique_ptr<const State3d> ConstPtr;

    typedef enum {
        POSITION_X = 0,
        POSITION_Y, // 1
        POSITION_Z, // 2

        QUATERNION_W, // 3
        QUATERNION_X, // 4
        QUATERNION_Y, // 5
        QUATERNION_Z, // 6

        LIN_VEL_X, // 7
        LIN_VEL_Y, // 8
        LIN_VEL_Z, // 9

        ANG_VEL_X, // 10
        ANG_VEL_Y, // 11
        ANG_VEL_Z, // 12

        LIN_ACC_X, // 13
        LIN_ACC_Y, // 14
        LIN_ACC_Z, // 15

        ACC_BIAS_X, // 16
        ACC_BIAS_Y, // 17
        ACC_BIAS_Z, // 18

        GYR_BIAS_X, // 19
        GYR_BIAS_Y, // 20
        GYR_BIAS_Z, // 21

        N_STATES // 22
    } StateIdx_e;

    typedef enum {
        POS_ERR_X = 0,
        POS_ERR_Y, // 1
        POS_ERR_Z, // 2

        ATT_ERR_X, // 3
        ATT_ERR_Y, // 4
        ATT_ERR_Z, // 5

        LIN_VEL_ERR_X, // 6
        LIN_VEL_ERR_Y, // 7
        LIN_VEL_ERR_Z, // 8

        ANG_VEL_ERR_X, // 9
        ANG_VEL_ERR_Y, // 10
        ANG_VEL_ERR_Z, // 11

        LIN_ACC_ERR_X, // 12
        LIN_ACC_ERR_Y, // 13
        LIN_ACC_ERR_Z, // 14

        ACC_BIAS_ERR_X, // 15
        ACC_BIAS_ERR_Y, // 16
        ACC_BIAS_ERR_Z, // 17

        GYR_BIAS_ERR_X, // 18
        GYR_BIAS_ERR_Y, // 19
        GYR_BIAS_ERR_Z, // 20

        N_ERROR_STATES // 21
    } ErrorStateIdx_e;

    typedef enum {
        ANG_VEL_NOISE_X = 0,
        ANG_VEL_NOISE_Y,
        ANG_VEL_NOISE_Z,

        LIN_ACC_NOISE_X,  // 3
        LIN_ACC_NOISE_Y,
        LIN_ACC_NOISE_Z,

        ACC_BIAS_NOISE_X,  // 6
        ACC_BIAS_NOISE_Y,
        ACC_BIAS_NOISE_Z,

        GYR_BIAS_NOISE_X,  // 9
        GYR_BIAS_NOISE_Y,
        GYR_BIAS_NOISE_Z,

        N_PROCESS_NOISES  // 12
    } ProcessNoiseIdx_e;

    State3d() : State() {
//        PRINT_DEBUG("State3d::State3d")
        x_ = VectorX::Zero(N_STATES);
        P_ = MatrixX::Zero(N_ERROR_STATES, N_ERROR_STATES);
    }  ///< Default constructor

    State3d(double _tStamp) : State(_tStamp) {
//        PRINT_DEBUG("State3d::State3d @ t=" << timestamp_)
        x_ = VectorX::Zero(N_STATES);
        P_ = MatrixX::Zero(N_ERROR_STATES, N_ERROR_STATES);
    }

    State3d(const VectorX& _x, const MatrixX& _P, const double& _tStamp) : State(_x, _P, _tStamp) {
//        PRINT_DEBUG("State3d::State3d @ t=" << timestamp_ << ", with x=" << x_.transpose())
        assert(x_.size() == N_STATES);
        assert(P_.rows() == N_ERROR_STATES && P_.cols() == N_ERROR_STATES);
    }

  virtual ~State3d() {
//  PRINT_DEBUG("State3d::~State3d")
  }  ///< Default destructor

  void checkStateVector(void) {
    assert(x_.size() == N_STATES);
    State::checkStateVector();
  }

  void checkCovarianceMatrix(void) {
    assert(P_.rows() == N_ERROR_STATES && P_.cols() == N_ERROR_STATES);
    State::checkCovarianceMatrix();
  }

  void checkStateAndCovariance(void) {
    assert(x_.size() == N_STATES);
    assert(P_.rows() == N_ERROR_STATES && P_.cols() == N_ERROR_STATES);
    State::checkStateAndCovariance();
  }
};  // class State3d

}  // namespace ekf
