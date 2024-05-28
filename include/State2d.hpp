///
/// \file State2d.hpp
/// \brief Base class for an EKF state in 2D
/// \author Lionel GENEVE
/// \date 17/11/2018
/// \version 1.0
///
#pragma once

// ---------- Headers ----------
// STD
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

class State2d : public State {
 public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    typedef std::unique_ptr<State2d> Ptr;
    typedef std::unique_ptr<const State2d> ConstPtr;

    typedef enum {
        POSITION_X = 0,
        POSITION_Y, // 1

        ATTITUDE_Z, // 2

        LIN_VEL_X, // 3
        LIN_VEL_Y, // 4

        ANG_VEL_Z, // 5

        LIN_ACC_X, // 6
        LIN_ACC_Y, // 7

        ACC_BIAS_X, // 8
        ACC_BIAS_Y, // 9

        GYR_BIAS_Z, // 10

        N_STATES // 11
    } StateIdx_e;

    typedef enum {
        POS_ERR_X = 0,
        POS_ERR_Y, // 1

        ATT_ERR_Z, // 2

        LIN_VEL_ERR_X, // 3
        LIN_VEL_ERR_Y, // 4

        ANG_VEL_ERR_Z, // 5

        LIN_ACC_ERR_X, // 6
        LIN_ACC_ERR_Y, // 7

        ACC_BIAS_ERR_X, // 8
        ACC_BIAS_ERR_Y, // 9

        GYR_BIAS_ERR_Z, // 10

        N_ERROR_STATES // 11
    } ErrorStateIdx_e;

    typedef enum {
      ANG_VEL_NOISE_Z = 0,

      LIN_ACC_NOISE_X,  // 1
      LIN_ACC_NOISE_Y,  // 2

      ACC_BIAS_NOISE_X,  // 3
      ACC_BIAS_NOISE_Y,  // 4

      GYR_BIAS_NOISE_Z,  // 5

      N_PROCESS_NOISES  // 6
    } ProcessNoiseIdx_e;

    State2d() : State() {
        //PRINT_DEBUG("State2d::State2d")
        x_ = VectorX::Zero(N_STATES);
        P_ = MatrixX::Zero(N_ERROR_STATES, N_ERROR_STATES);
    }  ///< Default constructor

    State2d(double _tStamp) : State(_tStamp) {
        //PRINT_DEBUG("State2d::State2d @ t=" << timestamp_)
        x_ = VectorX::Zero(N_STATES);
        P_ = MatrixX::Zero(N_ERROR_STATES, N_ERROR_STATES);
    }

    State2d(const VectorX& _x, const MatrixX& _P, const double& _tStamp) : State(_x, _P, _tStamp) {
        //PRINT_DEBUG("State2d::State2d @ t=" << timestamp_ << ", with x=" << x_.transpose())
        assert(x_.size() == N_STATES);
        assert(P_.rows() == N_ERROR_STATES && P_.cols() == N_ERROR_STATES);
    }

  virtual ~State2d() {
//    PRINT_DEBUG("State2d::~State2d")
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
};  // class State2d

}  // namespace ekf
