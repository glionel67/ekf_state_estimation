///
/// \file Ekf2d.cpp
/// \brief 2D state estimation using an Extended Kalman filter.
/// \author Lionel GENEVE
/// \date 17/11/2018
/// \version 1.0
///

// ---------- Headers ----------
// STD/C++ headers
#include <iostream>
#include <cmath>
#include <cassert>

// Project headers
#include "Ekf2d.hpp"

namespace ekf {

Ekf2d::Ekf2d(const FilterParameters& _params) : Ekf<State2d>(_params) {
    PRINT_DEBUG("Ekf2d::Ekf2d")

    I_ = MatrixX::Identity(State2d::N_ERROR_STATES, State2d::N_ERROR_STATES);
    F_ = MatrixX::Identity(State2d::N_STATES, State2d::N_STATES);
}

void Ekf2d::init(const VectorX& _x0, const MatrixX& _P0, double _tStamp) {
    PRINT_DEBUG("Ekf2d::init")

    if (_x0.size() != State2d::N_STATES) {
        throw std::invalid_argument("Ekf2d::init: state vector size error: _x0 != N_STATES");
    }

    if (_P0.rows() != State2d::N_ERROR_STATES && _P0.cols() != State2d::N_ERROR_STATES) {
        throw std::invalid_argument(
            "Ekf2d::init: state covariance matrix size error: _P0.rows() || _P0.cols() != N_ERROR_STATES");
    }

    Ekf<State2d>::init(_x0, _P0, _tStamp);
}  // init

void Ekf2d::predict(const State2d* _currState, State2d* _nextState) {
    //PRINT_DEBUG("Ekf2d::predict:\n-from: " << *_currState << "\n-to: " << *_nextState)

    double dt = _nextState->timestamp() - _currState->timestamp();  // Time interval [s]
    if (dt <= 0.) {
        throw std::runtime_error("Ekf2d::predict: dt <= 0.");
    }
    //PRINT_DEBUG("Ekf2d::predict: dt=" << dt << " s")

    const VectorX xCurr = _currState->x(); // Get current state vector
    double heading = xCurr(State2d::ATTITUDE_Z);
    double vx = xCurr(State2d::LIN_VEL_X);
    double vy = xCurr(State2d::LIN_VEL_Y);
    double ax = xCurr(State2d::LIN_ACC_X);
    double ay = xCurr(State2d::LIN_ACC_Y);
    double wz = xCurr(State2d::ANG_VEL_Z);

    // Predict next state vector
    Matrix2 R;
    R << std::cos(heading), -std::sin(heading),
        std::sin(heading), std::cos(heading);

    // Position: p_{k+1} = p_{k} + R_{k} * (v_{k} * dt + 0.5 * dt^2 * a_{k})
    _nextState->x().segment(State2d::POSITION_X, 2) = xCurr.segment(State2d::POSITION_X, 2)
        + R * (xCurr.segment(State2d::LIN_VEL_X, 2) * dt + 0.5 * xCurr.segment(State2d::LIN_ACC_X, 2) * dt * dt);

    // Attitude: q_{k+1} = q_{k} \otimes q({\omega_{k} * dt)
    _nextState->x()(State2d::ATTITUDE_Z) = piToPi(heading + wz * dt);  // Keep angle in (-pi, pi]

    // Linear velocity: v_{k+1} = v_{k} + a_{k} * dt
    _nextState->x().segment(State2d::LIN_VEL_X, 2) = xCurr.segment(State2d::LIN_VEL_X, 2) + xCurr.segment(State2d::LIN_ACC_X, 2) * dt;

    // Angular velocity: \omega_{k+1} = \omega_{k}
    _nextState->x()(State2d::ANG_VEL_Z) = wz;

    // Linear acceleration: a_{k+1} = a_{k}
    _nextState->x().segment(State2d::LIN_ACC_X, 2) = xCurr.segment(State2d::LIN_ACC_X, 2);

    // Accelero bias: ba_{k+1} = ba_{k}
    _nextState->x().segment(State2d::ACC_BIAS_X, 2) = xCurr.segment(State2d::ACC_BIAS_X, 2);

    // Gyro bias: bg_{k+1} = bg_{k}
    _nextState->x()(State2d::GYR_BIAS_Z) = xCurr(State2d::GYR_BIAS_Z);

    // State transition matrix

    // dp/dp
    //F_.block<2, 2>(State2d::POS_ERR_X, State2d::POS_ERR_X) = Matrix2::Identity();
    // dp/dt
    F_(State2d::POS_ERR_X, State2d::ATT_ERR_Z) = -(vx * std::sin(heading) + vy * std::cos(heading)) * dt 
        - 0.5 * (ax * std::sin(heading) + ay * std::cos(heading)) * dt * dt;
    F_(State2d::POS_ERR_Y, State2d::ATT_ERR_Z) = (vx * std::cos(heading) - vy * std::sin(heading)) * dt
        + 0.5 * (ax * std::cos(heading) - ay * std::sin(heading)) * dt * dt;
    // dp/dv
    F_.block<2, 2>(State2d::POS_ERR_X, State2d::LIN_VEL_ERR_X) = R * dt;
    // dp/da
    F_.block<2, 2>(State2d::POS_ERR_X, State2d::LIN_ACC_ERR_X) = 0.5 * R * dt * dt;

    // dt/dt
    //F_(State2d::ATT_ERR_Z, State2d::ATT_ERR_Z) = 1.0;
    // dt/dw
    F_(State2d::ATT_ERR_Z, State2d::ANG_VEL_ERR_Z) = dt;

    // dv/dv
    //F_(State2d::LIN_VEL_ERR_X, State2d::LIN_VEL_ERR_X) = 1.0;
    F_(State2d::LIN_VEL_ERR_X, State2d::LIN_VEL_ERR_Y) = wz * dt;
    F_(State2d::LIN_VEL_ERR_Y, State2d::LIN_VEL_ERR_X) = -wz * dt;
    //F_(State2d::LIN_VEL_ERR_Y, State2d::LIN_VEL_ERR_Y) = 1.0;
    // dv/dw
    F_(State2d::LIN_VEL_ERR_X, State2d::ANG_VEL_ERR_Z) = vy * dt;
    F_(State2d::LIN_VEL_ERR_Y, State2d::ANG_VEL_ERR_Z) = -vx * dt;
    // dv/da
    F_(State2d::LIN_VEL_ERR_X, State2d::LIN_ACC_ERR_X) = dt;
    F_(State2d::LIN_VEL_ERR_Y, State2d::LIN_ACC_ERR_Y) = dt;

    // dw/dw
    //F_(State2d::ANG_VEL_ERR_Z, State2d::ANG_VEL_ERR_Z) = 1.0;

    // da/da
    //F_.block<2, 2>(State2d::LIN_ACC_ERR_X, State2d::LIN_ACC_ERR_X) = Matrix2::Identity();

    // dba/dba
    //F_.block<2, 2>(State2d::ACC_BIAS_ERR_X, State2d::ACC_BIAS_ERR_X) = Matrix2::Identity();

    // dbg/dbg
    //F_(State2d::GYR_BIAS_ERR_Z, State2d::GYR_BIAS_ERR_Z) = 1.0;
    
    assert(!F_.array().isNaN().any());
    assert(F_.allFinite());

    //PRINT_DEBUG("Ekf2d::predict: F=" << F_)

    // Process noise covariance matrix
    // Don't add process noise when the robot is static ???
    //MatrixX Q;
    //MatrixX W;

    // Propagate covariance matrix: P = F * P * F.t + W * Q * W.t;
    //_nextState->P() = F_ * _currState->P() * F_.transpose() + W * Q * W.transpose();
    _nextState->P() = F_ * _currState->P() * F_.transpose();
    _nextState->P()(State2d::ANG_VEL_ERR_Z, State2d::ANG_VEL_ERR_Z) += std::pow(filterParameters_.angVelProcessNoiseStd, 2.0); // std::pow((1e-1 * M_PI / 180.) / 3., 2.);
    _nextState->P().block<2, 2>(State2d::LIN_ACC_ERR_X, State2d::LIN_ACC_ERR_X) += Matrix2::Identity() * std::pow(filterParameters_.linAccProcessNoiseStd, 2.0); // std::pow(1e-3 / 3., 2.);
    _nextState->P().block<2, 2>(State2d::ACC_BIAS_ERR_X, State2d::ACC_BIAS_ERR_X) += Matrix2::Identity() * std::pow(filterParameters_.accBiasNoiseStd, 2.0);
    _nextState->P()(State2d::GYR_BIAS_ERR_Z, State2d::GYR_BIAS_ERR_Z) += std::pow(filterParameters_.gyrBiasNoiseStd, 2.0);

    // Check state vector and covariance matrix
    _nextState->checkStateAndCovariance();
}  // predict

bool Ekf2d::correctAccelerometer(State2d* /*_state*/, const Measurement* /*_meas*/) {
  return true; 
}  // correctAccelerometer

bool Ekf2d::correctGyroscope(State2d* /*_state*/, const Measurement* /*_meas*/) {
  return true;
}  // correctGyroscope

bool Ekf2d::correctLinVel(State2d* _state, const Measurement* _meas) {
    //PRINT_DEBUG("Ekf2d::correctLinVel: for " << *_state << " with " << *_meas)

    const int measSize = 2;
    assert(_meas->z().size() == measSize);
    assert(_meas->R().rows() == measSize && _meas->R().cols() == measSize);

    VectorX x = _state->x();
    VectorX deltaX = VectorX::Zero(State2d::N_ERROR_STATES);

    // Measurement function Jacobian
    MatrixX H = MatrixX::Zero(measSize, State2d::N_ERROR_STATES);
    H.block<2, 2>(0, State2d::LIN_VEL_ERR_X) = Matrix2::Identity();

    // Measurement noise (co)variance
    MatrixX R = _meas->R();

    // Predict measurement from state vector
    VectorX zPred = VectorX::Zero(measSize);

    // Innovation vector
    VectorX nu = VectorX::Zero(measSize);

    MatrixX PHt;

    // Innovation matrix
    MatrixX S;
    MatrixX Sinv;

    // Kalman gain
    MatrixX K;

    for (size_t n = 0; n < filterParameters_.nIters; n++) {
        // Predict measurement from state vector
        zPred = x.segment(State2d::LIN_VEL_X, measSize);
        assert(zPred.allFinite());

        // Innovation vector
        nu = _meas->z() - zPred;

        // Compute measurement Jacobian matrix
        //H.block<2, 2>(0, State2d::LIN_VEL_ERR_X) = Matrix2::Identity();

        // Innovation matrix
        PHt = _state->P() * H.transpose();
        S = H * PHt + R;
        if (!S.allFinite()) {
            std::cerr << "Ekf2d::correctLinVel: S=\n" << S << " problem\n";
            return false;
        }

        // Mahalanobis outlier test
        Sinv = S.inverse();
        double maha = nu.transpose() * Sinv * nu;
        //assert(maha >= 0.);
        if (maha < 0.) {
            std::cerr << "Ekf2d::correctLinVel: maha=" << maha << " < 0\n";
            return false;
        }
        //maha = std::sqrt(maha);
        if (maha > _meas->mahaThresh()) {
            std::cerr << "Ekf2d::correctLinVel: maha=" << maha
                    << " > mahaThresh=" << _meas->mahaThresh() << std::endl;
            return false;
        }

        // Kalman gain
        K = PHt * Sinv;

        // State vector correction
        deltaX = K * nu;
        assert(deltaX.allFinite());
        addErrorState(x, deltaX); // x += deltaX;

        // Stop if no further improvement
        if (deltaX.norm() < 1e-6)
            break;
    }
    
    // State vector correction
    assert(x.allFinite());
    _state->x() = x;
    
    // Covariance matrix correction
    MatrixX IKH = I_ - K * H;
    _state->P() = IKH * _state->P() * IKH.transpose() + K * R * K.transpose();
    errorReset(_state->P(), deltaX);

    // Check state vector and covariance matrix
    _state->checkStateAndCovariance();
    
    return true;
}  // correctLinVel

bool Ekf2d::correctAngVel(State2d* _state, const Measurement* _meas) {
    //PRINT_DEBUG("Ekf2d::correctAngVel: for " << *_state << " with " << *_meas)

    assert(_meas->z().size() == 1);
    assert(_meas->R().rows() == 1 && _meas->R().cols() == 1);

    VectorX x = _state->x();
    VectorX deltaX = VectorX::Zero(State2d::N_ERROR_STATES);

    // Measurement function Jacobian
    MatrixX H = MatrixX::Zero(1, State2d::N_ERROR_STATES);
    H(0, State2d::ANG_VEL_ERR_Z) = 1.0;

    // Measurement noise (co)variance
    MatrixX R = _meas->R();

    // Predict measurement from state vector
    VectorX zPred = VectorX::Zero(1);

    // Innovation vector
    VectorX nu = VectorX::Zero(1);

    MatrixX PHt;

    // Innovation matrix
    MatrixX S;
    MatrixX Sinv;

    // Kalman gain
    MatrixX K;

    for (size_t n = 0; n < filterParameters_.nIters; n++) {
        // Predict measurement from state vector
        zPred = x.segment(State2d::ANG_VEL_Z, 1);
        assert(zPred.allFinite());

        // Innovation vector
        nu = _meas->z() - zPred;

        // Compute measurement Jacobian matrix
        //H(0, State2d::ANG_VEL_ERR_Z) = 1.0;

        // Innovation matrix
        PHt = _state->P() * H.transpose();
        S = H * PHt + R;
        if (!S.allFinite()) {
            std::cerr << "Ekf2d::correctAngVel: S=\n" << S << " problem\n";
            return false;
        }

        // Mahalanobis outlier test
        Sinv = S.inverse();
        double maha = nu.transpose() * Sinv * nu;
        //assert(maha >= 0.);
        if (maha < 0.) {
            std::cerr << "Ekf2d::correctAngVel: maha=" << maha << " < 0\n";
            return false;
        }
        //maha = std::sqrt(maha);
        if (maha > _meas->mahaThresh()) {
            std::cerr << "Ekf2d::correctAngVel: maha=" << maha
                    << " > mahaThresh=" << _meas->mahaThresh() << std::endl;
            return false;
        }

        // Kalman gain
        K = PHt * Sinv;

        // State vector correction
        deltaX = K * nu;
        assert(deltaX.allFinite());
        addErrorState(x, deltaX); // x += deltaX;

        // Stop if no further improvement
        if (deltaX.norm() < 1e-6)
            break;
    }
    
    // State vector correction
    assert(x.allFinite());
    _state->x() = x;
    
    // Covariance matrix correction
    MatrixX IKH = I_ - K * H;
    _state->P() = IKH * _state->P() * IKH.transpose() + K * R * K.transpose();
    errorReset(_state->P(), deltaX);

    // Check state vector and covariance matrix
    _state->checkStateAndCovariance();
    
    return true;
}  // correctLinVel

bool Ekf2d::correctLinAngVel(State2d* _state, const Measurement* _meas) {
    //PRINT_DEBUG("Ekf2d::correctLinAngVel: for " << *_state << " with " << *_meas)

    const int measSize = 3;
    assert(_meas->z().size() == measSize);
    assert(_meas->R().rows() == measSize && _meas->R().cols() == measSize);

    VectorX x = _state->x();  // State vector
    VectorX deltaX = VectorX::Zero(State2d::N_ERROR_STATES);  // Error state vector

    // Measurement function Jacobian
    MatrixX H = MatrixX::Zero(measSize, State2d::N_ERROR_STATES);
    H.block<2, 2>(0, State2d::LIN_VEL_ERR_X) = Matrix2::Identity();
    H(2, State2d::ANG_VEL_ERR_Z) = 1.0;

    // Measurement noise (co)variance
    MatrixX R = _meas->R();

    // Predict measurement from state vector
    VectorX zPred = VectorX::Zero(measSize);

    // Innovation vector
    VectorX nu = VectorX::Zero(measSize);

    MatrixX PHt;

    // Innovation matrix
    MatrixX S;
    MatrixX Sinv;

    // Kalman gain
    MatrixX K;

    for (size_t n = 0; n < filterParameters_.nIters; n++) {
        // Predict
        zPred.segment(0, 2) = x.segment(State2d::LIN_VEL_X, 2);
        zPred(2) = x(State2d::ANG_VEL_Z);
        assert(zPred.allFinite());
        assert(!zPred.array().isNaN().any());
        //PRINT_DEBUG("Ekf2d::correctLinAngVel: z=" << _meas->z().transpose() << ", zPred=" << zPred.transpose())

        // Innovation vector
        nu = _meas->z() - zPred;
        //PRINT_DEBUG("Ekf2d::correctLinAngVel: nu=" << nu.transpose())

        // Compute measurement Jacobian matrix
        //H.block<2, 2>(0, State2d::LIN_VEL_ERR_X) = Matrix2::Identity();
        //H(2, State2d::ANG_VEL_ERR_Z) = 1.0;
        //PRINT_DEBUG("Ekf2d::correctLinAngVel: size(H)=(" << H.rows() << ", " << H.cols() 
        //        << "), size(P)=(" << _state->P().rows() << ", " << _state->P().cols() << ")")

        // Innovation matrix
        PHt = _state->P() * H.transpose();
        S = H * PHt + R;
        if (!S.allFinite()) {
            std::cerr << "Ekf2d::correctLinAngVel: S=\n" << S << " problem\n";
            return false;
        }

        // Mahalanobis outlier test
        Sinv = S.inverse();
        double maha = nu.transpose() * Sinv * nu;
        //assert(maha >= 0.);
        if (maha < 0.) {
            std::cerr << "Ekf2d::correctLinAngVel: maha=" << maha << " < 0\n";
            return false;
        }
        //maha = std::sqrt(maha);
        if (maha > _meas->mahaThresh()) {
            std::cerr << "Ekf2d::correctLinAngVel: maha=" << maha << " > thresh=" << _meas->mahaThresh() << std::endl;
            // We still continue applying the measurement but with much higher covariance (This is to avoid divergence of the solution).
            if (maha <= 3.0 * _meas->mahaThresh()) {
                S = H * PHt + R * (maha / _meas->mahaThresh());
                Sinv = S.inverse();
            } else {
                return false;
            }
        }

        // Kalman gain
        K = PHt * Sinv;

        // State vector correction
        deltaX = K * nu;
        assert(deltaX.allFinite());
        assert(!deltaX.array().isNaN().any());
        addErrorState(x, deltaX); // x += deltaX;

        // Stop if no further improvement
        if (deltaX.norm() < 1e-6)
            break;
    }
    
    // State vector correction
    assert(x.allFinite());
    assert(!x.array().isNaN().any());
    _state->x() = x;
    
    // Covariance matrix correction
    MatrixX IKH = I_ - K * H;
    _state->P() = IKH * _state->P() * IKH.transpose() + K * R * K.transpose();
    errorReset(_state->P(), deltaX);

    // Check state vector and covariance matrix
    _state->checkStateAndCovariance();
    
    return true;
}  // correctLinAngVel

bool Ekf2d::correctLinAcc(State2d* _state, const Measurement* _meas) {
    //PRINT_DEBUG("Ekf2d::correctLinAcc: for " << *_state << " with " << *_meas)

    const int measSize = 2;
    assert(_meas->z().size() == measSize);
    assert(_meas->R().rows() == measSize && _meas->R().cols() == measSize);

    VectorX x = _state->x();
    VectorX deltaX = VectorX::Zero(State2d::N_ERROR_STATES);

    // Measurement function Jacobian
    MatrixX H = MatrixX::Zero(measSize, State2d::N_ERROR_STATES);
    H.block<2, 2>(0, State2d::LIN_ACC_ERR_X) = Matrix2::Identity();

    // Measurement noise (co)variance
    MatrixX R = _meas->R();

    // Predict measurement from state vector
    VectorX zPred = VectorX::Zero(measSize);

    // Innovation vector
    VectorX nu = VectorX::Zero(measSize);

    MatrixX PHt;

    // Innovation matrix
    MatrixX S;
    MatrixX Sinv;

    // Kalman gain
    MatrixX K;

    for (size_t n = 0; n < filterParameters_.nIters; n++) {
        // Predict measurement from state vector
        zPred = x.segment(State2d::LIN_ACC_X, measSize);
        assert(zPred.allFinite());
        assert(!zPred.array().isNaN().any());

        // Innovation vector
        nu = _meas->z() - zPred;

        // Compute measurement Jacobian matrix
        //H.block<2, 2>(0, State2d::LIN_ACC_ERR_X) = Matrix2::Identity();

        // Innovation matrix
        PHt = _state->P() * H.transpose();
        S = H * PHt + R;
        if (!S.allFinite()) {
            std::cerr << "Ekf2d::correctLinAcc: S=\n" << S << " problem\n";
            return false;
        }

        // Mahalanobis outlier test
        Sinv = S.inverse();
        double maha = nu.transpose() * Sinv * nu;
        //assert(maha >= 0.);
        if (maha < 0.) {
            std::cerr << "Ekf2d::correctLinAcc: maha=" << maha << " < 0\n";
            return false;
        }
        //maha = std::sqrt(maha);
        if (maha > _meas->mahaThresh()) {
            std::cerr << "Ekf2d::correctLinAcc: maha=" << maha << " > thresh=" << _meas->mahaThresh() << std::endl;
            // We still continue applying the measurement but with much higher covariance (This is to avoid divergence of the solution).
            if (maha <= 3.0 * _meas->mahaThresh()) {
                S = H * PHt + R * (maha / _meas->mahaThresh());
                Sinv = S.inverse();
            } else {
                return false;
            }
        }

        // Kalman gain
        K = PHt * Sinv;

        // State vector correction
        deltaX = K * nu;
        assert(deltaX.allFinite());
        assert(!deltaX.array().isNaN().any());
        addErrorState(x, deltaX); // x += deltaX;

        // Stop if no further improvement
        if (deltaX.norm() < 1e-6)
            break;
    }
    
    // State vector correction
    assert(x.allFinite());
    assert(!x.array().isNaN().any());
    _state->x() = x;
    
    // Covariance matrix correction
    MatrixX IKH = I_ - K * H;
    _state->P() = IKH * _state->P() * IKH.transpose() + K * R * K.transpose();
    errorReset(_state->P(), deltaX);

    // Check state vector and covariance matrix
    _state->checkStateAndCovariance();
    
    return true;
}  // correctLinAcc

bool Ekf2d::correctPosition(State2d* _state, const Measurement* _meas) {
    //PRINT_DEBUG("Ekf2d::correctPosition: for " << *_state << " with " << *_meas)

    const int measSize = 2;
    assert(_meas->z().size() == measSize);
    assert(_meas->R().rows() == measSize && _meas->R().cols() == measSize);

    VectorX x = _state->x();  // State vector
    VectorX deltaX = VectorX::Zero(State2d::N_ERROR_STATES);  // Error-state vector

    // Measurement function Jacobian
    MatrixX H = MatrixX::Zero(measSize, State2d::N_ERROR_STATES);
    H.block<2, 2>(0, State2d::POS_ERR_X) = Matrix2::Identity();

    // Measurement noise (co)variance
    MatrixX R = _meas->R();

    // Predict measurement from state vector
    VectorX zPred = VectorX::Zero(measSize);

    // Innovation vector
    VectorX nu = VectorX::Zero(measSize);

    MatrixX PHt;

    // Innovation matrix
    MatrixX S;
    MatrixX Sinv;

    // Kalman gain
    MatrixX K;

    for (size_t n = 0; n < filterParameters_.nIters; n++) {
        // Predict
        zPred = x.segment(State2d::POSITION_X, measSize);
        assert(zPred.allFinite());
        assert(!zPred.array().isNaN().any());

        // Innovation vector
        nu = _meas->z() - zPred;

        // Compute measurement Jacobian matrix
        //H.block<2, 2>(0, State2d::POS_ERR_X) = Matrix2::Identity();

        // Innovation matrix
        PHt = _state->P() * H.transpose();
        S = H * PHt + R;
        if (!S.allFinite()) {
            std::cerr << "Ekf2d::correctPosition: S=\n" << S << " problem\n";
            return false;
        }

        // Mahalanobis outlier test
        Sinv = S.inverse();
        double maha = nu.transpose() * Sinv * nu;
        //assert(maha >= 0.);
        if (maha < 0.) {
            std::cerr << "Ekf2d::correctPosition: maha=" << maha << " < 0\n";
            return false;
        }
        //maha = std::sqrt(maha);
        if (maha > _meas->mahaThresh()) {
            std::cerr << "Ekf2d::correctPosition: maha=" << maha << " > thresh=" << _meas->mahaThresh() << std::endl;
            // We still continue applying the measurement but with much higher covariance (This is to avoid divergence of the solution).
            if (maha <= 3.0 * _meas->mahaThresh()) {
                S = H * PHt + R * (maha / _meas->mahaThresh());
                Sinv = S.inverse();
            } else {
                return false;
            }
        }

        // Kalman gain
        K = PHt * Sinv;

        // State vector correction
        deltaX = K * nu;
        assert(deltaX.allFinite());
        assert(!deltaX.array().isNaN().any());
        addErrorState(x, deltaX); // x += deltaX;

        // Stop if no further improvement
        if (deltaX.norm() < 1e-6)
            break;
    }
    
    // State vector correction
    assert(x.allFinite());
    assert(!x.array().isNaN().any());
    _state->x() = x;
    
    // Covariance matrix correction
    MatrixX IKH = I_ - K * H;
    _state->P() = IKH * _state->P() * IKH.transpose() + K * R * K.transpose();
    errorReset(_state->P(), deltaX);

    // Check state vector and covariance matrix
    _state->checkStateAndCovariance();
    
    return true;
}  // correctPosition

bool Ekf2d::correctAttitude(State2d* _state, const Measurement* _meas) {
    //PRINT_DEBUG("Ekf2d::correctAttitude: for " << *_state << " with " << *_meas)

    VectorX x = _state->x();
    VectorX deltaX = VectorX::Zero(State2d::N_ERROR_STATES);

    // Measurement function Jacobian
    MatrixX H = MatrixX::Zero(1, State2d::N_ERROR_STATES);
    H(0, State2d::ATT_ERR_Z) = 1.0;

    // Measurement noise (co)variance
    MatrixX R = _meas->R();

    // Predict measurement from state vector
    VectorX zPred = VectorX::Zero(1);

    // Innovation vector
    VectorX nu = VectorX::Zero(1);

    MatrixX PHt;

    // Innovation matrix
    MatrixX S;
    MatrixX Sinv;

    // Kalman gain
    MatrixX K;

    for (size_t n = 0; n < filterParameters_.nIters; n++) {
        // Predict
        zPred(0) = x(State2d::ATTITUDE_Z);
        assert(zPred.allFinite());

        // Innovation vector
        nu = _meas->z() - zPred;

        // Compute measurement Jacobian matrix
        //H(0, State2d::ATT_ERR_Z) = 1.0;

        // Innovation matrix
        PHt = _state->P() * H.transpose();
        S = H * PHt + R;
        if (!S.allFinite()) {
            std::cerr << "Ekf2d::correctAttitude: S=\n" << S << " problem\n";
            return false;
        }

        // Mahalanobis outlier test
        Sinv = S.inverse();
        double maha = nu.transpose() * Sinv * nu;
        //assert(maha >= 0.);
        if (maha < 0.) {
            std::cerr << "Ekf2d::correctAttitude: maha=" << maha << " < 0\n";
            return false;
        }
        //maha = std::sqrt(maha);
        if (maha > _meas->mahaThresh()) {
            std::cerr << "Ekf2d::correctAttitude: maha=" << maha
                    << " > mahaThresh=" << _meas->mahaThresh() << std::endl;
            return false;
        }

        // Kalman gain
        K = PHt * Sinv;

        // State vector correction
        deltaX = K * nu;
        assert(deltaX.allFinite());
        addErrorState(x, deltaX); // x += deltaX;

        // Stop if no further improvement
        if (deltaX.norm() < 1e-6)
            break;
    }
    
    // State vector correction
    assert(x.allFinite());
    _state->x() = x;
    
    // Covariance matrix correction
    MatrixX IKH = I_ - K * H;
    _state->P() = IKH * _state->P() * IKH.transpose() + K * R * K.transpose();
    errorReset(_state->P(), deltaX);

    // Check state vector and covariance matrix
    _state->checkStateAndCovariance();
    
    return true;
}

bool Ekf2d::correctAccGyr(State2d* _state, const Measurement* _meas) {
    //PRINT_DEBUG("Ekf2d::correctAccGyr: for " << *_state << " with " << *_meas)

    const int measSize = 3;
    assert(_state->x().size() == State2d::N_STATES);
    assert(_state->P().rows() == State2d::N_ERROR_STATES && _state->P().cols() == State2d::N_ERROR_STATES);
    assert(_meas->z().size() == measSize);
    assert(_meas->R().rows() == measSize && _meas->R().cols() == measSize);

    VectorX x = _state->x();  // State vector
    VectorX deltaX = VectorX::Zero(State2d::N_ERROR_STATES);  // Error state vector

    // Measurement function Jacobian
    MatrixX H = MatrixX::Zero(measSize, State2d::N_ERROR_STATES);
    H.block<2, 2>(0, State2d::LIN_ACC_ERR_X) = Matrix2::Identity();
    H.block<2, 2>(0, State2d::ACC_BIAS_ERR_X) = Matrix2::Identity();
    H(2, State2d::ANG_VEL_ERR_Z) = 1.0;
    H(2, State2d::GYR_BIAS_ERR_Z) = 1.0;

    // Measurement noise (co)variance
    MatrixX R = _meas->R();

    // Predict measurement from state vector
    VectorX zPred = VectorX::Zero(measSize);

    // Innovation vector
    VectorX nu = VectorX::Zero(measSize);

    MatrixX PHt;

    // Innovation matrix
    MatrixX S;
    MatrixX Sinv;

    // Kalman gain
    MatrixX K;

    for (size_t n = 0; n < filterParameters_.nIters; n++) {
        // Predict measurement from current state
        zPred.segment(0, 2) = x.segment(State2d::LIN_ACC_X, 2) + x.segment(State2d::ACC_BIAS_X, 2);  // Acc
        zPred(2) = x(State2d::ANG_VEL_Z) + x(State2d::GYR_BIAS_Z);  // Gyr
        assert(zPred.allFinite());
        assert(!zPred.array().isNaN().any());
        //PRINT_DEBUG("Ekf2d::correctAccGyr: z=" << _meas->z().transpose() << ", zPred=" << zPred.transpose());

        // Innovation vector
        nu = _meas->z() - zPred;
        //PRINT_DEBUG("Ekf2d::correctAccGyr: nu=" << nu.transpose())

        // Compute measurement Jacobian matrix
        //H.block<2, 2>(0, State2d::LIN_ACC_ERR_X) = Matrix2::Identity();
        //H.block<2, 2>(0, State2d::ACC_BIAS_ERR_X) = Matrix2::Identity();
        //H(2, State2d::ANG_VEL_ERR_Z) = 1.0;
        //H(2, State2d::GYR_BIAS_ERR_Z) = 1.0;
        //PRINT_DEBUG("Ekf2d::correctAccGyr: size(H)=(" << H.rows() << ", " << H.cols() 
        //        << "), size(P)=(" << _state->P().rows() << ", " << _state->P().cols() << ")");

        // Innovation matrix
        PHt = _state->P() * H.transpose();
        S = H * PHt + R;
        //PRINT_DEBUG("Ekf2d::correctAccGyr: size(S)=(" << S.rows() << ", " << S.cols() << ")")
        if (!S.allFinite()) {
            std::cerr << "Ekf2d::correctAccGyr: S=\n" << S << " problem\n";
            return false;
        }

        // Mahalanobis outlier test
        Sinv = S.inverse();
        double maha = nu.transpose() * Sinv * nu;
        //PRINT_DEBUG("Ekf2d::correctAccGyr: maha=" << maha) 
        //assert(maha >= 0.);
        if (maha < 0.) {
            std::cerr << "Ekf2d::correctAccGyr: maha=" << maha << " < 0\n";
            return false;
        }
        //maha = std::sqrt(maha);
        if (maha > _meas->mahaThresh()) {
            std::cerr << "Ekf2d::correctAccGyr: maha=" << maha << " > thresh=" << _meas->mahaThresh() << std::endl;
            // We still continue applying the measurement but with much higher covariance (This is to avoid divergence of the solution).
            if (maha <= 3.0 * _meas->mahaThresh()) {
                S = H * PHt + R * (maha / _meas->mahaThresh());
                Sinv = S.inverse();
            } else {
                return false;
            }
        }

        // Kalman gain
        K = PHt * Sinv;
        //PRINT_DEBUG("Ekf2d::correctAccGyr: size(K)=(" << K.rows() << ", " << K.cols() << ")")

        // State vector correction
        deltaX = K * nu;
        assert(deltaX.allFinite());
        assert(!deltaX.array().isNaN().any());
        addErrorState(x, deltaX); // x += deltaX;

        // Stop if no further improvement
        if (deltaX.norm() < 1e-6)
            break;
    }
    
    // State vector correction
    assert(x.allFinite());
    assert(!x.array().isNaN().any());
    _state->x() = x;
    
    // Covariance matrix correction
    MatrixX IKH = I_ - K * H;
    _state->P() = IKH * _state->P() * IKH.transpose() + K * R * K.transpose();
    errorReset(_state->P(), deltaX);  // Reset error-state

    // Check state vector and covariance matrix
    _state->checkStateAndCovariance();
    
    return true;
}  // correctAccGyr

void Ekf2d::addErrorState(VectorX& _x, VectorX& _dx) {
    //PRINT_DEBUG("Ekf2d::addErrorState");
    // 2D position
    _x.segment(State2d::POSITION_X, 2) += _dx.segment(State2d::POS_ERR_X, 2);
    // 2D attitude
    _x(State2d::ATTITUDE_Z) += _dx(State2d::ATT_ERR_Z);
    // 2D linear velocity
    _x.segment(State2d::LIN_VEL_X, 2) += _dx.segment(State2d::LIN_VEL_ERR_X, 2);
    // 2D angular velocity
    _x(State2d::ANG_VEL_Z) += _dx(State2d::ANG_VEL_ERR_Z);
    // 2D linear acceleration
    _x.segment(State2d::LIN_ACC_X, 2) += _dx.segment(State2d::LIN_ACC_ERR_X, 2);
    // Accelero bias
    _x.segment(State2d::ACC_BIAS_X, 2) += _dx.segment(State2d::ACC_BIAS_ERR_X, 2);
    // Gyro bias
    _x(State2d::GYR_BIAS_Z) += _dx(State2d::GYR_BIAS_ERR_Z);
}  // addErrorState

void Ekf2d::errorReset(MatrixX& _P, VectorX& _dx) {
    //PRINT_DEBUG("Ekf2d::errorReset");
    assert(_P.rows() == State2d::N_ERROR_STATES && _P.cols() == State2d::N_ERROR_STATES);
    assert(_dx.size() == State2d::N_ERROR_STATES);

    MatrixX G = MatrixX::Identity(State2d::N_ERROR_STATES, State2d::N_ERROR_STATES);
    _P = G * _P * G.transpose();
}  // errorReset

bool Ekf2d::checkSanity(void) {
    //PRINT_DEBUG("Ekf2d::checkSanity");
    // Check filter state vector
    stateSanityMask_ = 0;
    VectorX lastState = states_.back().x();
    if (lastState(State2d::POSITION_X) > filterParameters_.maxPosition(0)
    || lastState(State2d::POSITION_X) < filterParameters_.minPosition(0)
    || lastState(State2d::POSITION_Y) > filterParameters_.maxPosition(1)
    || lastState(State2d::POSITION_Y) < filterParameters_.minPosition(1)) {
        stateSanityMask_ |= 0b1;
    }
    if (lastState(State2d::ATTITUDE_Z) < filterParameters_.minAttitude(2)
    || lastState(State2d::ATTITUDE_Z) > filterParameters_.maxAttitude(2)) {
        stateSanityMask_ |= 0b10;
    }
    if (lastState(State2d::LIN_VEL_X) > filterParameters_.maxLinVel(0)
    || lastState(State2d::LIN_VEL_X) < filterParameters_.minLinVel(0)
    || lastState(State2d::LIN_VEL_Y) > filterParameters_.maxLinVel(1)
    || lastState(State2d::LIN_VEL_Y) < filterParameters_.minLinVel(1)) {
        stateSanityMask_ |= 0b100;
    }
    if (lastState(State2d::ANG_VEL_Z) > filterParameters_.maxAngVel(2)
    || lastState(State2d::ANG_VEL_Z) < filterParameters_.minAngVel(2)) {
        stateSanityMask_ |= 0b1000;
    }
    if (lastState(State2d::LIN_ACC_X) > filterParameters_.maxLinAcc(0)
    || lastState(State2d::LIN_ACC_X) < filterParameters_.minLinAcc(0)
    || lastState(State2d::LIN_ACC_Y) > filterParameters_.maxLinAcc(1)
    || lastState(State2d::LIN_ACC_Y) < filterParameters_.minLinAcc(1)) {
        stateSanityMask_ |= 0b10000;
    }
    if (lastState(State2d::ACC_BIAS_X) > filterParameters_.maxAccBias(0)
    || lastState(State2d::ACC_BIAS_X) < filterParameters_.minAccBias(0)
    || lastState(State2d::ACC_BIAS_Y) > filterParameters_.maxAccBias(1)
    || lastState(State2d::ACC_BIAS_Y) < filterParameters_.minAccBias(1)) {
        stateSanityMask_ |= 0b100000;
    }
    if (lastState(State2d::GYR_BIAS_Z) > filterParameters_.maxGyrBias(2)
    || lastState(State2d::GYR_BIAS_Z) < filterParameters_.minGyrBias(2)) {
        stateSanityMask_ |= 0b1000000;
    }

    // Check filter covariance matrix
    covSanityMask_ = 0;
    MatrixX lastCov = states_.back().P();
    if (lastCov(State2d::POS_ERR_X, State2d::POS_ERR_X) > std::pow(filterParameters_.maxPositionStd(0), 2.0)
    || lastCov(State2d::POS_ERR_Y, State2d::POS_ERR_Y) > std::pow(filterParameters_.maxPositionStd(1), 2.0)) {
        covSanityMask_ |= 0b1;
    }
    if (lastCov(State2d::ATT_ERR_Z, State2d::ATT_ERR_Z) > std::pow(filterParameters_.maxAttitudeStd(2), 2.0)) {
        covSanityMask_ |= 0b10;
    }
    if (lastCov(State2d::LIN_VEL_ERR_X, State2d::LIN_VEL_ERR_X) > std::pow(filterParameters_.maxLinVelStd(0), 2.0)
    || lastCov(State2d::LIN_VEL_ERR_Y, State2d::LIN_VEL_ERR_Y) > std::pow(filterParameters_.maxLinVelStd(1), 2.0)) {
        covSanityMask_ |= 0b100;
    }
    if (lastCov(State2d::ANG_VEL_ERR_Z, State2d::ANG_VEL_ERR_Z) > std::pow(filterParameters_.maxAngVelStd(2), 2.0)) {
        covSanityMask_ |= 0b1000;
    }
    if (lastCov(State2d::LIN_ACC_ERR_X, State2d::LIN_ACC_ERR_X) > std::pow(filterParameters_.maxLinAccStd(0), 2.0)
    || lastCov(State2d::LIN_ACC_ERR_Y, State2d::LIN_ACC_ERR_Y) > std::pow(filterParameters_.maxLinAccStd(1), 2.0)) {
        covSanityMask_ |= 0b10000;
    }
    if (lastCov(State2d::ACC_BIAS_ERR_X, State2d::ACC_BIAS_ERR_X) > std::pow(filterParameters_.maxAccBiasStd(0), 2.0)
    || lastCov(State2d::ACC_BIAS_ERR_Y, State2d::ACC_BIAS_ERR_Y) > std::pow(filterParameters_.maxAccBiasStd(1), 2.0)) {
        covSanityMask_ |= 0b100000;
    }
    if (lastCov(State2d::GYR_BIAS_ERR_Z, State2d::GYR_BIAS_ERR_Z) > std::pow(filterParameters_.maxGyrBiasStd(2), 2.0)) {
        covSanityMask_ |= 0b1000000;
    }

    isSane_ = (stateSanityMask_ > 0) || (covSanityMask_ > 0) ? false : true;

    return isSane_;
}  // checkSanity

} // namespace ekf
