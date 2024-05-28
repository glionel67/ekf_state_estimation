///
/// \file Ekf3d.cpp
/// \brief 3D state estimation using an Extended Kalman filter.
/// Useful link on IMU: https://www.vectornav.com/resources/inertial-navigation-primer/math-fundamentals/math-coning#:~:text=The%20coning%20integral%20provides%20a,amount%20of%20time%2C%20%CE%94t.
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
#include "Ekf3d.hpp"


namespace ekf {

Ekf3d::Ekf3d(const FilterParameters& _params) : Ekf<State3d>(_params) {
    PRINT_DEBUG("Ekf3d::Ekf3d")

    I_ = MatrixX::Identity(State3d::N_ERROR_STATES, State3d::N_ERROR_STATES);
    F_ = MatrixX::Identity(State3d::N_STATES, State3d::N_STATES);
}

void Ekf3d::init(const VectorX& _x0, const MatrixX& _P0, double _tStamp) {
    PRINT_DEBUG("Ekf3d::init")

    if (_x0.size() != State3d::N_STATES) {
        throw std::invalid_argument("Ekf3d::init: state vector size error: _x0 != N_STATES");
    }

    if (_P0.rows() != State3d::N_ERROR_STATES && _P0.cols() != State3d::N_ERROR_STATES) {
        throw std::invalid_argument(
            "Ekf3d::init: state covariance matrix size error: _P0.rows() || _P0.cols() != N_ERROR_STATES");
    }

    Ekf<State3d>::init(_x0, _P0, _tStamp);
}  // init

void Ekf3d::predict(const State3d* _currState, State3d* _nextState) {
    //PRINT_DEBUG("Ekf3d::predict:\n-from: " << *_currState << "\n-to: " << *_nextState)

    double dt = _nextState->timestamp() - _currState->timestamp();  // Time interval [s]
    //assert(dt > 0.)
    if (dt <= 0.) {
        throw std::runtime_error("Ekf3d::predict: dt <= 0.");
    }
    //PRINT_DEBUG("Ekf3d::predict: dt=" << dt << " s")

    const VectorX xCurr = _currState->x(); // Get current state vector

    // Predict next state vector
    Quaternion qCurr(xCurr(State3d::QUATERNION_W), xCurr(State3d::QUATERNION_X),
                    xCurr(State3d::QUATERNION_Y), xCurr(State3d::QUATERNION_Z));
    Matrix3 R = qCurr.toRotationMatrix();

    // Position: p_{k+1} = p_{k} + R_{k} * (v_{k} * dt + 0.5 * dt^2 * a_{k})
    _nextState->x().segment(State3d::POSITION_X, 3) = xCurr.segment(State3d::POSITION_X, 3)
        + R * (xCurr.segment(State3d::LIN_VEL_X, 3) * dt + 0.5 * xCurr.segment(State3d::LIN_ACC_X, 3) * dt * dt);

    // Attitude: q_{k+1} = q_{k} \otimes q({\omega_{k} * dt)
    Vector3 deltaAngle = xCurr.segment(State3d::ANG_VEL_X, 3) * dt;
    Quaternion qRot = deltaAngleToQuaternion(deltaAngle);
    Quaternion qNext = qCurr * qRot;
    qNext.normalize();  // Normalize quaternion
    _nextState->x()(State3d::QUATERNION_W) = qNext.w();
    _nextState->x()(State3d::QUATERNION_X) = qNext.x();
    _nextState->x()(State3d::QUATERNION_Y) = qNext.y();
    _nextState->x()(State3d::QUATERNION_Z) = qNext.z();

    // Linear velocity: v_{k+1} = v_{k} + a_{k} * dt
    _nextState->x().segment(State3d::LIN_VEL_X, 3) = xCurr.segment(State3d::LIN_VEL_X, 3) + xCurr.segment(State3d::LIN_ACC_X, 3) * dt;

    // Angular velocity: \omega_{k+1} = \omega_{k}
    _nextState->x().segment(State3d::ANG_VEL_X, 3) = xCurr.segment(State3d::ANG_VEL_X, 3);

    // Linear acceleration: a_{k+1} = a_{k}
    _nextState->x().segment(State3d::LIN_ACC_X, 3) = xCurr.segment(State3d::LIN_ACC_X, 3);

    // Accelero bias: ba_{k+1} = ba_{k}
    _nextState->x().segment(State3d::ACC_BIAS_X, 3) = xCurr.segment(State3d::ACC_BIAS_X, 3);

    // Gyro bias: bg_{k+1} = bg_{k}
    _nextState->x().segment(State3d::GYR_BIAS_X, 3) = xCurr.segment(State3d::GYR_BIAS_X, 3);

    // State transition matrix

    // dp/dp
    //F.block<3, 3>(State3d::POS_ERR_X, State3d::POS_ERR_X) = Matrix3::Identity();
    // dp/dt
    Matrix3 Vx;
    Vx << 0., -xCurr(State3d::LIN_VEL_Z), xCurr(State3d::LIN_VEL_Y),
        xCurr(State3d::LIN_VEL_Z), 0., -xCurr(State3d::LIN_VEL_X),
        -xCurr(State3d::LIN_VEL_Y), xCurr(State3d::LIN_VEL_X), 0.;
    F_.block<3, 3>(State3d::POS_ERR_X, State3d::ATT_ERR_X) = -R * Vx * dt;
    // dp/dv
    F_.block<3, 3>(State3d::POS_ERR_X, State3d::LIN_VEL_ERR_X) = R * dt;

    // dt/dt
    F_.block<3, 3>(State3d::ATT_ERR_X, State3d::ATT_ERR_X) = qRot.toRotationMatrix();
    // dt/dw
    F_.block<3, 3>(State3d::ATT_ERR_X, State3d::ANG_VEL_ERR_X) = dt * Matrix3::Identity();

    // dv/dv
    //F_.block<3, 3>(State3d::LIN_VEL_ERR_X, State3d::LIN_VEL_ERR_X) = Matrix3::Identity();
    // dv/da
    F_.block<3, 3>(State3d::LIN_VEL_ERR_X, State3d::LIN_ACC_ERR_X) = dt * Matrix3::Identity();

    // dw/dw
    //F_.block<3, 3>(State3d::ANG_VEL_ERR_X, State3d::ANG_VEL_ERR_X) = Matrix3::Identity();

    // da/da
    //F_.block<3, 3>(State3d::LIN_ACC_ERR_X, State3d::LIN_ACC_ERR_X) = Matrix3::Identity();

    // dba/dba
    //F_.block<3, 3>(State3d::ACC_BIAS_ERR_X, State3d::ACC_BIAS_ERR_X) = Matrix3::Identity();

    // dbg/dbg
    //F_.block<3, 3>(State3d::GYR_BIAS_ERR_X, State3d::GYR_BIAS_ERR_X) = Matrix3::Identity();
    
    assert(!F_.array().isNaN().any());
    assert(F_.allFinite());

    //PRINT_DEBUG("Ekf3d::predict: F=" << F_)

    // Process noise covariance matrix
    // TODO: Don't add process noise when the robot is static ???
    // Version 1
    //MatrixX Q = MatrixX::Zero(State3d::N_ERROR_STATES, State3d::N_ERROR_STATES);
    //Q.block<3, 3>(State3d::POS_ERR_X, State3d::POS_ERR_X) = Matrix3::Identity() * std::pow(filterParameters_.positionProcessNoiseStd, 2.0);  // Position --- std::pow(1e-4 / 3., 2.)
    //Q.block<3, 3>(State3d::ATT_ERR_X, State3d::ATT_ERR_X) = Matrix3::Identity() * std::pow(filterParameters_.attitudeProcessNoiseStd, 2.0);  // Attitude --- std::pow((.1 * M_PI / 180.) / 3., 2.)
    //Q.block<3, 3>(State3d::LIN_VEL_ERR_X, State3d::LIN_VEL_ERR_X) = Matrix3::Identity() * std::pow(filterParameters_.linVelProcessNoiseStd, 2.0);  // Linear velocity --- std::pow(1e-3 / 3., 2.)
    //Q.block<3, 3>(State3d::ANG_VEL_ERR_X, State3d::ANG_VEL_ERR_X) = Matrix3::Identity() * std::pow(filterParameters_.angVelProcessNoiseStd, 2.0);  // Angular velocity --- std::pow((1e-1 * M_PI / 180.) / 3., 2.)
    //Q.block<3, 3>(State3d::LIN_ACC_ERR_X, State3d::LIN_ACC_ERR_X) = Matrix3::Identity() * std::pow(filterParameters_.linAccProcessNoiseStd, 2.0); // Linear acceleration --- std::pow(1e-3 / 3., 2.)
    //Q.block<3, 3>(State3d::ACC_BIAS_ERR_X, State3d::ACC_BIAS_ERR_X) = Matrix3::Identity() * std::pow(filterParameters_.accBiasNoiseStd, 2.0);  // Accelero bias
    //Q.block<3, 3>(State3d::GYR_BIAS_ERR_X, State3d::GYR_BIAS_ERR_X) = Matrix3::Identity() * std::pow(filterParameters_.gyrBiasNoiseStd, 2.0); // Gyro bias
    // Version 2
    //MatrixX Q = MatrixX::Zero(State3d::N_PROCESS_NOISES, State3d::N_PROCESS_NOISES);
    //Q.block<3, 3>(State3d::ANG_VEL_NOISE_X, State3d::ANG_VEL_NOISE_X) = Matrix3::Identity() * std::pow(filterParameters_.angVelProcessNoiseStd, 2.0); // std::pow((1e-1 * M_PI / 180.) / 3., 2.);
    //Q.block<3, 3>(State3d::LIN_ACC_NOISE_X, State3d::LIN_ACC_NOISE_X) = Matrix3::Identity() * std::pow(filterParameters_.linAccProcessNoiseStd, 2.0); // std::pow(1e-3 / 3., 2.);
    //Q.block<3, 3>(State3d::ACC_BIAS_NOISE_X, State3d::ACC_BIAS_NOISE_X) = Matrix3::Identity() * std::pow(filterParameters_.accBiasNoiseStd, 2.0);
    //Q.block<3, 3>(State3d::GYR_BIAS_NOISE_X, State3d::GYR_BIAS_NOISE_X) = Matrix3::Identity() * std::pow(filterParameters_.gyrBiasNoiseStd, 2.0);
    //assert(!Q.array().isNaN().any());
    //assert(Q.allFinite());
    //PRINT_DEBUG("Ekf3d::predict: Q=" << Q)

    //MatrixX W = MatrixX::Zero(State3d::N_ERROR_STATES, State3d::N_PROCESS_NOISES);
    //W.block<3, 3>(State3d::ANG_VEL_ERR_X, State3d::ANG_VEL_NOISE_X) = Matrix3::Identity();
    //W.block<3, 3>(State3d::LIN_ACC_ERR_X, State3d::LIN_ACC_NOISE_X) = Matrix3::Identity();
    //W.block<3, 3>(State3d::ACC_BIAS_ERR_X, State3d::ACC_BIAS_NOISE_X) = Matrix3::Identity();
    //W.block<3, 3>(State3d::GYR_BIAS_ERR_X, State3d::GYR_BIAS_NOISE_X) = Matrix3::Identity();
    //assert(!W.array().isNaN().any());
    //assert(W.allFinite());
    //PRINT_DEBUG("Ekf3d::predict: W=" << W)

    // Propagate covariance matrix: P = F * P * F.t + W * Q * W.t;
    //_nextState->P() = F_ * _currState->P() * F_.transpose() + W * Q * W.transpose();
    _nextState->P() = F_ * _currState->P() * F_.transpose();
    _nextState->P().block<3, 3>(State3d::ANG_VEL_ERR_X, State3d::ANG_VEL_ERR_X) += Matrix3::Identity() * std::pow(filterParameters_.angVelProcessNoiseStd, 2.0); // std::pow((1e-1 * M_PI / 180.) / 3., 2.);
    _nextState->P().block<3, 3>(State3d::LIN_ACC_ERR_X, State3d::LIN_ACC_ERR_X) += Matrix3::Identity() * std::pow(filterParameters_.linAccProcessNoiseStd, 2.0); // std::pow(1e-3 / 3., 2.);
    _nextState->P().block<3, 3>(State3d::ACC_BIAS_ERR_X, State3d::ACC_BIAS_ERR_X) += Matrix3::Identity() * std::pow(filterParameters_.accBiasNoiseStd, 2.0);
    _nextState->P().block<3, 3>(State3d::GYR_BIAS_ERR_X, State3d::GYR_BIAS_ERR_X) += Matrix3::Identity() * std::pow(filterParameters_.gyrBiasNoiseStd, 2.0);

    // Check state vector and covariance matrix
    _nextState->checkStateAndCovariance();
}  // predict

bool Ekf3d::correctLinVel(State3d* _state, const Measurement* _meas) {
    //PRINT_DEBUG("Ekf3d::correctLinVel: for " << *_state << " with " << *_meas)

    assert(_meas->z().size() == 3);
    assert(_meas->R().rows() == 3 && _meas->R().cols() == 3);

    VectorX x = _state->x();
    VectorX deltaX = VectorX::Zero(State3d::N_ERROR_STATES);

    // Measurement function Jacobian
    MatrixX H = MatrixX::Zero(3, State3d::N_ERROR_STATES);
    H.block<3, 3>(0, State3d::LIN_VEL_ERR_X) = Matrix3::Identity();

    // Measurement noise (co)variance
    MatrixX R = _meas->R();

    // Predict measurement from state vector
    VectorX zPred = VectorX::Zero(3);

    // Innovation vector
    VectorX nu = VectorX::Zero(3);

    MatrixX PHt;

    // Innovation matrix
    MatrixX S;
    MatrixX Sinv;

    // Kalman gain
    MatrixX K;

    for (size_t n = 0; n < filterParameters_.nIters; n++) {
        // Predict measurement from state vector
        zPred = x.segment(State3d::LIN_VEL_X, 3);
        assert(zPred.allFinite());

        // Innovation vector
        nu = _meas->z() - zPred;

        // Compute measurement Jacobian matrix
        //H.block<3, 3>(0, LIN_VEL_ERR_X) = Matrix3::Identity();

        // Innovation matrix
        PHt = _state->P() * H.transpose();
        S = H * PHt + R;
        if (!S.allFinite()) {
            std::cerr << "Ekf3d::correctLinVel: S=\n" << S << " problem\n";
            return false;
        }

        // Mahalanobis outlier test
        Sinv = S.inverse();
        double maha = nu.transpose() * Sinv * nu;
        //assert(maha >= 0.);
        if (maha < 0.) {
            std::cerr << "Ekf3d::correctLinVel: maha=" << maha << " < 0\n";
            return false;
        }
        //maha = std::sqrt(maha);
        if (maha > _meas->mahaThresh()) {
            std::cerr << "Ekf3d::correctLinVel: maha=" << maha
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

bool Ekf3d::correctAngVel(State3d* _state, const Measurement* _meas) {
    //PRINT_DEBUG("Ekf3d::correctAngVel: for " << *_state << " with " << *_meas)

    assert(_meas->z().size() == 3);
    assert(_meas->R().rows() == 3 && _meas->R().cols() == 3);

    VectorX x = _state->x();
    VectorX deltaX = VectorX::Zero(State3d::N_ERROR_STATES);

    // Measurement function Jacobian
    MatrixX H = MatrixX::Zero(3, State3d::N_ERROR_STATES);
    H.block<3, 3>(0, State3d::ANG_VEL_ERR_X) = Matrix3::Identity();

    // Measurement noise (co)variance
    MatrixX R = _meas->R();

    // Predict measurement from state vector
    VectorX zPred = VectorX::Zero(3);

    // Innovation vector
    VectorX nu = VectorX::Zero(3);

    MatrixX PHt;

    // Innovation matrix
    MatrixX S;
    MatrixX Sinv;

    // Kalman gain
    MatrixX K;

    for (size_t n = 0; n < filterParameters_.nIters; n++) {
        // Predict measurement from state vector
        zPred = x.segment(State3d::ANG_VEL_X, 3);
        assert(zPred.allFinite());

        // Innovation vector
        nu = _meas->z() - zPred;

        // Compute measurement Jacobian matrix
        //H.block<3, 3>(0, ANG_VEL_ERR_X) = Matrix3::Identity();

        // Innovation matrix
        PHt = _state->P() * H.transpose();
        S = H * PHt + R;
        if (!S.allFinite()) {
            std::cerr << "Ekf3d::correctAngVel: S=\n" << S << " problem\n";
            return false;
        }

        // Mahalanobis outlier test
        Sinv = S.inverse();
        double maha = nu.transpose() * Sinv * nu;
        //assert(maha >= 0.);
        if (maha < 0.) {
            std::cerr << "Ekf3d::correctAngVel: maha=" << maha << " < 0\n";
            return false;
        }
        //maha = std::sqrt(maha);
        if (maha > _meas->mahaThresh()) {
            std::cerr << "Ekf3d::correctAngVel: maha=" << maha
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

bool Ekf3d::correctLinAngVel(State3d* _state, const Measurement* _meas) {
    //PRINT_DEBUG("Ekf3d::correctLinAngVel: for " << *_state << " with " << *_meas)

    assert(_meas->z().size() == 6);
    assert(_meas->R().rows() == 6 && _meas->R().cols() == 6);

    VectorX x = _state->x();  // State vector
    VectorX deltaX = VectorX::Zero(State3d::N_ERROR_STATES);  // Error state vector

    // Measurement function Jacobian
    MatrixX H = MatrixX::Zero(6, State3d::N_ERROR_STATES);
    H.block<3, 3>(0, State3d::LIN_VEL_ERR_X) = Matrix3::Identity();
    H.block<3, 3>(3, State3d::ANG_VEL_ERR_X) = Matrix3::Identity();

    // Measurement noise (co)variance
    MatrixX R = _meas->R();

    // Predict measurement from state vector
    VectorX zPred = VectorX::Zero(6);

    // Innovation vector
    VectorX nu = VectorX::Zero(6);

    MatrixX PHt;

    // Innovation matrix
    MatrixX S;
    MatrixX Sinv;

    // Kalman gain
    MatrixX K;

    for (size_t n = 0; n < filterParameters_.nIters; n++) {
        // Predict
        zPred.segment(0, 3) = x.segment(State3d::LIN_VEL_X, 3);
        zPred.segment(3, 3) = x.segment(State3d::ANG_VEL_X, 3);
        assert(zPred.allFinite());
        assert(!zPred.array().isNaN().any());
        //PRINT_DEBUG("Ekf3d::correctLinAngVel: z=" << _meas->z().transpose() << ", zPred=" << zPred.transpose())

        // Innovation vector
        nu = _meas->z() - zPred;
        //PRINT_DEBUG("Ekf3d::correctLinAngVel: nu=" << nu.transpose())

        // Compute measurement Jacobian matrix
        //H.block<3, 3>(0, State3d::LIN_VEL_ERR_X) = Matrix3::Identity();
        //H.block<3, 3>(3, State3d::ANG_VEL_ERR_X) = Matrix3::Identity();
        //PRINT_DEBUG("Ekf3d::correctLinAngVel: size(H)=(" << H.rows() << ", " << H.cols() 
        //        << "), size(P)=(" << _state->P().rows() << ", " << _state->P().cols() << ")")

        // Innovation matrix
        PHt = _state->P() * H.transpose();
        S = H * PHt + R;
        if (!S.allFinite()) {
            std::cerr << "Ekf3d::correctLinAngVel: S=\n" << S << " problem\n";
            return false;
        }

        // Mahalanobis outlier test
        Sinv = S.inverse();
        double maha = nu.transpose() * Sinv * nu;
        //assert(maha >= 0.);
        if (maha < 0.) {
            std::cerr << "Ekf3d::correctLinAngVel: maha=" << maha << " < 0\n";
            return false;
        }
        //maha = std::sqrt(maha);
        if (maha > _meas->mahaThresh()) {
            std::cerr << "Ekf3d::correctLinAngVel: maha=" << maha << " > thresh=" << _meas->mahaThresh() << std::endl;
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

bool Ekf3d::correctLinAcc(State3d* _state, const Measurement* _meas) {
    //PRINT_DEBUG("Ekf3d::correctLinAcc: for " << *_state << " with " << *_meas)

    assert(_meas->z().size() == 3);
    assert(_meas->R().rows() == 3 && _meas->R().cols() == 3);

    VectorX x = _state->x();
    VectorX deltaX = VectorX::Zero(State3d::N_ERROR_STATES);

    // Measurement function Jacobian
    MatrixX H = MatrixX::Zero(3, State3d::N_ERROR_STATES);
    H.block<3, 3>(0, State3d::LIN_ACC_ERR_X) = Matrix3::Identity();

    // Measurement noise (co)variance
    MatrixX R = _meas->R();

    // Predict measurement from state vector
    VectorX zPred = VectorX::Zero(3);

    // Innovation vector
    VectorX nu = VectorX::Zero(3);

    MatrixX PHt;

    // Innovation matrix
    MatrixX S;
    MatrixX Sinv;

    // Kalman gain
    MatrixX K;

    for (size_t n = 0; n < filterParameters_.nIters; n++) {
        // Predict measurement from state vector
        zPred = x.segment(State3d::LIN_ACC_X, 3);
        assert(zPred.allFinite());
        assert(!zPred.array().isNaN().any());

        // Innovation vector
        nu = _meas->z() - zPred;

        // Compute measurement Jacobian matrix
        //H.block<3, 3>(0, State3d::LIN_ACC_ERR_X) = Matrix3::Identity();

        // Innovation matrix
        PHt = _state->P() * H.transpose();
        S = H * PHt + R;
        if (!S.allFinite()) {
            std::cerr << "Ekf3d::correctLinAcc: S=\n" << S << " problem\n";
            return false;
        }

        // Mahalanobis outlier test
        Sinv = S.inverse();
        double maha = nu.transpose() * Sinv * nu;
        //assert(maha >= 0.);
        if (maha < 0.) {
            std::cerr << "Ekf3d::correctLinAcc: maha=" << maha << " < 0\n";
            return false;
        }
        //maha = std::sqrt(maha);
        if (maha > _meas->mahaThresh()) {
            std::cerr << "Ekf3d::correctLinAcc: maha=" << maha << " > thresh=" << _meas->mahaThresh() << std::endl;
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

bool Ekf3d::correctPosition(State3d* _state, const Measurement* _meas) {
    //PRINT_DEBUG("Ekf3d::correctPosition: for " << *_state << " with " << *_meas)

    assert(_meas->z().size() == 3);
    assert(_meas->R().rows() == 3 && _meas->R().cols() == 3);

    VectorX x = _state->x();  // State vector
    VectorX deltaX = VectorX::Zero(State3d::N_ERROR_STATES);  // Error-state vector

    // Measurement function Jacobian
    MatrixX H = MatrixX::Zero(3, State3d::N_ERROR_STATES);
    H.block<3, 3>(0, State3d::POS_ERR_X) = Matrix3::Identity();

    // Measurement noise (co)variance
    MatrixX R = _meas->R();

    // Predict measurement from state vector
    VectorX zPred = VectorX::Zero(3);

    // Innovation vector
    VectorX nu = VectorX::Zero(3);

    MatrixX PHt;

    // Innovation matrix
    MatrixX S;
    MatrixX Sinv;

    // Kalman gain
    MatrixX K;

    for (size_t n = 0; n < filterParameters_.nIters; n++) {
        // Predict
        zPred = x.segment(State3d::POSITION_X, 3);
        assert(zPred.allFinite());
        assert(!zPred.array().isNaN().any());

        // Innovation vector
        nu = _meas->z() - zPred;

        // Compute measurement Jacobian matrix
        //H.block<3, 3>(0, State3d::POS_ERR_X) = Matrix3::Identity();

        // Innovation matrix
        PHt = _state->P() * H.transpose();
        S = H * PHt + R;
        if (!S.allFinite()) {
            std::cerr << "Ekf3d::correctPosition: S=\n" << S << " problem\n";
            return false;
        }

        // Mahalanobis outlier test
        Sinv = S.inverse();
        double maha = nu.transpose() * Sinv * nu;
        //assert(maha >= 0.);
        if (maha < 0.) {
            std::cerr << "Ekf3d::correctPosition: maha=" << maha << " < 0\n";
            return false;
        }
        //maha = std::sqrt(maha);
        if (maha > _meas->mahaThresh()) {
            std::cerr << "Ekf3d::correctPosition: maha=" << maha << " > thresh=" << _meas->mahaThresh() << std::endl;
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

bool Ekf3d::correctAttitude(State3d* _state, const Measurement* _meas) {
    //PRINT_DEBUG("Ekf3d::correctAttitude: for " << *_state << " with " << *_meas)

    VectorX x = _state->x();
    VectorX deltaX = VectorX::Zero(State3d::N_ERROR_STATES);

    // Measurement function Jacobian
    MatrixX H = MatrixX::Zero(3, State3d::N_ERROR_STATES);
    H.block<3, 3>(0, State3d::ATT_ERR_X) = Matrix3::Identity(); // TODO: replace with the correct formula

    // Measurement noise (co)variance
    MatrixX R = _meas->R();

    // Predict measurement from state vector
    VectorX zPred = VectorX::Zero(4);

    // Innovation vector
    VectorX nu = VectorX::Zero(4);

    MatrixX PHt;

    // Innovation matrix
    MatrixX S;
    MatrixX Sinv;

    // Kalman gain
    MatrixX K;

    for (size_t n = 0; n < filterParameters_.nIters; n++) {
        // Predict
        zPred = x.segment(State3d::QUATERNION_W, 4);
        assert(zPred.allFinite());

        // Innovation vector
        nu = _meas->z() - zPred;

        // Compute measurement Jacobian matrix
        //H.block<3, 3>(0, State3d::ATT_ERR_X) = Matrix3::Identity(); // TODO: replace with the correct formula

        // Innovation matrix
        PHt = _state->P() * H.transpose();
        S = H * PHt + R;
        if (!S.allFinite()) {
            std::cerr << "Ekf3d::correctAttitude: S=\n" << S << " problem\n";
            return false;
        }

        // Mahalanobis outlier test
        Sinv = S.inverse();
        double maha = nu.transpose() * Sinv * nu;
        //assert(maha >= 0.);
        if (maha < 0.) {
            std::cerr << "Ekf3d::correctAttitude: maha=" << maha << " < 0\n";
            return false;
        }
        //maha = std::sqrt(maha);
        if (maha > _meas->mahaThresh()) {
            std::cerr << "Ekf3d::correctAttitude: maha=" << maha
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

bool Ekf3d::correctAccelerometer(State3d* _state, const Measurement* _meas) {
    //PRINT_DEBUG("Ekf3d::correctAccelero: for " << *_state << " with " << *_meas)

    assert(_meas->z().size() == 3);
    assert(_meas->R().rows() == 3 && _meas->R().cols() == 3);

    VectorX x = _state->x();
    VectorX deltaX = VectorX::Zero(State3d::N_ERROR_STATES);

    // Measurement function Jacobian
    MatrixX H = MatrixX::Zero(3, State3d::N_ERROR_STATES);
    H.block<3, 3>(0, State3d::LIN_ACC_ERR_X) = Matrix3::Identity();
    H.block<3, 3>(0, State3d::ACC_BIAS_ERR_X) = Matrix3::Identity();

    // Measurement noise (co)variance
    MatrixX R = _meas->R();

    // Predict measurement from state vector
    VectorX zPred = VectorX::Zero(3);

    // Innovation vector
    VectorX nu = VectorX::Zero(3);

    MatrixX PHt;

    // Innovation matrix
    MatrixX S;
    MatrixX Sinv;

    // Kalman gain
    MatrixX K;

    for (size_t n = 0; n < filterParameters_.nIters; n++) {
        // Predict measurement from current state
        Quaternion q(x(State3d::QUATERNION_W), x(State3d::QUATERNION_X), 
                    x(State3d::QUATERNION_Y), x(State3d::QUATERNION_Z));
        Matrix3 Rot = q.toRotationMatrix();
        Vector3 Rtg = Rot.transpose() * filterParameters_.gravityVector;
        zPred = x.segment(State3d::LIN_ACC_X, 3) - Rtg + x.segment(State3d::ACC_BIAS_X, 3);
        assert(zPred.allFinite());

        // Innovation vector
        nu = _meas->z() - zPred;

        // Compute measurement Jacobian matrix
        Matrix3 Rtgx;
        Rtgx << 0.0, -Rtg(2), Rtg(1),
                Rtg(2), 0.0, -Rtg(0),
                -Rtg(1), Rtg(0), 0.0;
        H.block<3, 3>(0, State3d::ATT_ERR_X) = Rtgx;
        //H.block<3, 3>(0, State3d::LIN_ACC_ERR_X) = Matrix3::Identity();
        //H.block<3, 3>(0, State3d::ACC_BIAS_ERR_X) = Matrix3::Identity();

        // Innovation matrix
        PHt = _state->P() * H.transpose();
        S = H * PHt + R;
        if (!S.allFinite()) {
            std::cerr << "Ekf3d::correctAccelero: S=\n" << S << " problem\n";
            return false;
        }

        // Mahalanobis outlier test
        Sinv = S.inverse();
        double maha = nu.transpose() * Sinv * nu;
        //assert(maha >= 0.);
        if (maha < 0.) {
            std::cerr << "Ekf3d::correctAccelero: maha=" << maha << " < 0\n";
            return false;
        }
        //maha = std::sqrt(maha);
        if (maha > _meas->mahaThresh()) {
            std::cerr << "Ekf3d::correctAccelero: maha=" << maha
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
}  // correctAccelerometer

bool Ekf3d::correctGyroscope(State3d* _state, const Measurement* _meas) {
    //PRINT_DEBUG("Ekf3d::correctGyro: for " << *_state << " with " << *_meas)

    assert(_meas->z().size() == 3);
    assert(_meas->R().rows() == 3 && _meas->R().cols() == 3);

    VectorX x = _state->x();
    VectorX deltaX = VectorX::Zero(State3d::N_ERROR_STATES);

    // Measurement function Jacobian
    MatrixX H = MatrixX::Zero(3, State3d::N_ERROR_STATES);
    H.block<3, 3>(0, State3d::ANG_VEL_ERR_X) = Matrix3::Identity();
    H.block<3, 3>(0, State3d::GYR_BIAS_ERR_X) = Matrix3::Identity();

    // Measurement noise (co)variance
    MatrixX R = _meas->R();

    // Predict measurement from state vector
    VectorX zPred = VectorX::Zero(3);

    // Innovation vector
    VectorX nu = VectorX::Zero(3);

    MatrixX PHt;

    // Innovation matrix
    MatrixX S;
    MatrixX Sinv;

    // Kalman gain
    MatrixX K;

    for (size_t n = 0; n < filterParameters_.nIters; n++) {
        // Predict measurement from current state
        zPred = x.segment(State3d::ANG_VEL_X, 3) + x.segment(State3d::GYR_BIAS_X, 3);
        assert(zPred.allFinite());

        // Innovation vector
        nu = _meas->z() - zPred;

        // Compute measurement Jacobian matrix
        //H.block<3, 3>(0, State3d::ANG_VEL_ERR_X) = Matrix3::Identity();
        //H.block<3, 3>(0, State3d::GYR_BIAS_ERR_X) = Matrix3::Identity();

        // Innovation matrix
        PHt = _state->P() * H.transpose();
        S = H * PHt + R;
        if (!S.allFinite()) {
            std::cerr << "Ekf3d::correctGyro: S=\n" << S << " problem\n";
            return false;
        }

        // Mahalanobis outlier test
        Sinv = S.inverse();
        double maha = nu.transpose() * Sinv * nu;
        //assert(maha >= 0.);
        if (maha < 0.) {
            std::cerr << "Ekf3d::correctGyro: maha=" << maha << " < 0\n";
            return false;
        }
        //maha = std::sqrt(maha);
        if (maha > _meas->mahaThresh()) {
            std::cerr << "Ekf3d::correctGyro: maha=" << maha
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
}  // correctGyroscope

bool Ekf3d::correctAccGyr(State3d* _state, const Measurement* _meas) {
    //PRINT_DEBUG("Ekf3d::correctAccGyr: for " << *_state << " with " << *_meas)

    assert(_state->x().size() == State3d::N_STATES);
    assert(_state->P().rows() == State3d::N_ERROR_STATES && _state->P().cols() == State3d::N_ERROR_STATES);
    assert(_meas->z().size() == 6);
    assert(_meas->R().rows() == 6 && _meas->R().cols() == 6);

    VectorX x = _state->x();  // State vector
    VectorX deltaX = VectorX::Zero(State3d::N_ERROR_STATES);  // Error state vector

    // Measurement function Jacobian
    MatrixX H = MatrixX::Zero(6, State3d::N_ERROR_STATES);
    H.block<3, 3>(0, State3d::LIN_ACC_ERR_X) = Matrix3::Identity();
    H.block<3, 3>(0, State3d::ACC_BIAS_ERR_X) = Matrix3::Identity();
    H.block<3, 3>(3, State3d::ANG_VEL_ERR_X) = Matrix3::Identity();
    H.block<3, 3>(3, State3d::GYR_BIAS_ERR_X) = Matrix3::Identity();
    Matrix3 Rtgx;

    // Measurement noise (co)variance
    MatrixX R = _meas->R();

    // Predict measurement from state vector
    VectorX zPred = VectorX::Zero(6);

    // Innovation vector
    VectorX nu = VectorX::Zero(6);

    MatrixX PHt;

    // Innovation matrix
    MatrixX S;
    MatrixX Sinv;

    // Kalman gain
    MatrixX K;

    for (size_t n = 0; n < filterParameters_.nIters; n++) {
        // Predict measurement from current state
        Quaternion q(x(State3d::QUATERNION_W), x(State3d::QUATERNION_X),
                    x(State3d::QUATERNION_Y), x(State3d::QUATERNION_Z));
        Matrix3 Rot = q.toRotationMatrix();
        Vector3 Rtg = Rot.transpose() * filterParameters_.gravityVector;
        zPred.segment(0, 3) = x.segment(State3d::LIN_ACC_X, 3) + Rtg + x.segment(State3d::ACC_BIAS_X, 3);  // Acc
        zPred.segment(3, 3) = x.segment(State3d::ANG_VEL_X, 3) + x.segment(State3d::GYR_BIAS_X, 3);  // Gyr
        assert(zPred.allFinite());
        assert(!zPred.array().isNaN().any());
        //PRINT_DEBUG("Ekf3d::correctAccGyr: z=" << _meas->z().transpose() << ", zPred=" << zPred.transpose())

        // Innovation vector
        nu = _meas->z() - zPred;
        //PRINT_DEBUG("Ekf3d::correctAccGyr: nu=" << nu.transpose())

        // Compute measurement Jacobian matrix
        Rtgx << 0.0, -Rtg(2), Rtg(1),
                Rtg(2), 0.0, -Rtg(0),
                -Rtg(1), Rtg(0), 0.0;
        H.block<3, 3>(0, State3d::ATT_ERR_X) = Rtgx;
        //H.block<3, 3>(0, State3d::LIN_ACC_ERR_X) = Matrix3::Identity();
        //H.block<3, 3>(0, State3d::ACC_BIAS_ERR_X) = Matrix3::Identity();
        //H.block<3, 3>(3, State3d::ANG_VEL_ERR_X) = Matrix3::Identity();
        //H.block<3, 3>(3, State3d::GYR_BIAS_ERR_X) = Matrix3::Identity();
        //PRINT_DEBUG("Ekf3d::correctAccGyr: size(H)=(" << H.rows() << ", " << H.cols() 
        //        << "), size(P)=(" << _state->P().rows() << ", " << _state->P().cols() << ")")

        // Innovation matrix
        PHt = _state->P() * H.transpose();
        S = H * PHt + R;
        //PRINT_DEBUG("Ekf3d::correctAccGyr: size(S)=(" << S.rows() << ", " << S.cols() << ")")
        if (!S.allFinite()) {
            std::cerr << "Ekf3d::correctAccGyr: S=\n" << S << " problem\n";
            return false;
        }

        // Mahalanobis outlier test
        Sinv = S.inverse();
        double maha = nu.transpose() * Sinv * nu;
        //PRINT_DEBUG("Ekf3d::correctAccGyr: maha=" << maha) 
        //assert(maha >= 0.);
        if (maha < 0.) {
            std::cerr << "Ekf3d::correctAccGyr: maha=" << maha << " < 0\n";
            return false;
        }
        //maha = std::sqrt(maha);
        if (maha > _meas->mahaThresh()) {
            std::cerr << "Ekf3d::correctAccGyr: maha=" << maha << " > thresh=" << _meas->mahaThresh() << std::endl;
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
        //PRINT_DEBUG("Ekf3d::correctAccGyr: size(K)=(" << K.rows() << ", " << K.cols() << ")")

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

void Ekf3d::addErrorState(VectorX& _x, VectorX& _dx) {
    //PRINT_DEBUG("Ekf3d::addErrorState")
    // 3D position
    _x.segment(State3d::POSITION_X, 3) += _dx.segment(State3d::POS_ERR_X, 3);
    // 3D attitude
    Quaternion qCurr(_x(State3d::QUATERNION_W), _x(State3d::QUATERNION_X),
                    _x(State3d::QUATERNION_Y), _x(State3d::QUATERNION_Z));
    Vector3 deltaAngle = _dx.segment(State3d::ATT_ERR_X, 3);
    Quaternion qRot = deltaAngleToQuaternion(deltaAngle);
    Quaternion qNew = qCurr * qRot;
    qNew.normalize();  // Normalize quaternion
    _x(State3d::QUATERNION_W) = qNew.w();
    _x(State3d::QUATERNION_X) = qNew.x();
    _x(State3d::QUATERNION_Y) = qNew.y();
    _x(State3d::QUATERNION_Z) = qNew.z();
    // 3D linear velocity
    _x.segment(State3d::LIN_VEL_X, 3) += _dx.segment(State3d::LIN_VEL_ERR_X, 3);
    // 3D angular velocity
    _x.segment(State3d::ANG_VEL_X, 3) += _dx.segment(State3d::ANG_VEL_ERR_X, 3);
    // 3D linear acceleration
    _x.segment(State3d::LIN_ACC_X, 3) += _dx.segment(State3d::LIN_ACC_ERR_X, 3);
    // Accelero bias
    _x.segment(State3d::ACC_BIAS_X, 3) += _dx.segment(State3d::ACC_BIAS_ERR_X, 3);
    // Gyro bias
    _x.segment(State3d::GYR_BIAS_X, 3) += _dx.segment(State3d::GYR_BIAS_ERR_X, 3);
}  // addErrorState

void Ekf3d::errorReset(MatrixX& _P, VectorX& _dx) {
    assert(_P.rows() == State3d::N_ERROR_STATES && _P.cols() == State3d::N_ERROR_STATES);
    assert(_dx.size() == State3d::N_ERROR_STATES);

    // No rotation/attitude error, don't need to reset covariance matrix
    Vector3 attError = _dx.segment(State3d::ATT_ERR_X, 3);
    if (attError.norm() < 1e-9) {
        return;
    }

    // G_{att} = I_3 - [0.5 * \delta_{\theta}]_{\times}
    Matrix3 Gatt;
    Gatt << 0.0, -_dx(State3d::ATT_ERR_Z), _dx(State3d::ATT_ERR_Y),
        _dx(State3d::ATT_ERR_Z), 0.0, -_dx(State3d::ATT_ERR_X),
        -_dx(State3d::ATT_ERR_Y), _dx(State3d::ATT_ERR_X), 0.0;
    Gatt = Matrix3::Identity() - 0.5 * Gatt;

    MatrixX G = MatrixX::Identity(State3d::N_ERROR_STATES, State3d::N_ERROR_STATES);
    G.block<3, 3>(State3d::ATT_ERR_X, State3d::ATT_ERR_X) = Gatt;

    _P = G * _P * G.transpose();
}  // errorReset

Quaternion Ekf3d::deltaAngleToQuaternion(const Vector3& _deltaAngle) {
    double rotAngle = 0.5 * _deltaAngle.norm();  // Rotation angle
    Vector3 rotAxis = _deltaAngle;  // Rotation axis vector
    if (rotAngle >= 1e-5)  // Avoid small rotation problem
        rotAxis = _deltaAngle.normalized();

    double ca = cos(rotAngle);
    double sa = sin(rotAngle);
    
    return Quaternion(ca, rotAxis(0) * sa, rotAxis(1) * sa, rotAxis(2) * sa);
}  // deltaAngleToQuaternion

bool Ekf3d::checkSanity(void) {
    // Check filter state vector
    stateSanityMask_ = 0;  // Reset mask
    VectorX lastState = states_.back().x();

    // Check 3D position
    if (lastState(State3d::POSITION_X) > filterParameters_.maxPosition(0)
    || lastState(State3d::POSITION_X) < filterParameters_.minPosition(0)
    || lastState(State3d::POSITION_Y) > filterParameters_.maxPosition(1)
    || lastState(State3d::POSITION_Y) < filterParameters_.minPosition(1)
    || lastState(State3d::POSITION_Z) > filterParameters_.maxPosition(2)
    || lastState(State3d::POSITION_Z) < filterParameters_.minPosition(2)) {
        stateSanityMask_ |= 0b1;
    }

    // Check 3D attitude
    // Convert quaternion to roll/pitch/yaw Euler angles
    Quaternion q(lastState(State3d::QUATERNION_W), lastState(State3d::QUATERNION_X),
                lastState(State3d::QUATERNION_Y), lastState(State3d::QUATERNION_Z));
    double roll = 0.0;
    double pitch = 0.0;
    double yaw = 0.0;
    quaternionToRollPitchYaw(q, &roll, &pitch, &yaw);
    if (roll > filterParameters_.maxAttitude(0)
    || roll < filterParameters_.minAttitude(0)
    || pitch > filterParameters_.maxAttitude(1)
    || pitch < filterParameters_.minAttitude(1)
    || yaw > filterParameters_.maxAttitude(2)
    || yaw < filterParameters_.minAttitude(2)) {
      stateSanityMask_ |= 0b10;
    }

    // Check 3D linear velocity
    if (lastState(State3d::LIN_VEL_X) > filterParameters_.maxLinVel(0)
    || lastState(State3d::LIN_VEL_X) < filterParameters_.minLinVel(0)
    || lastState(State3d::LIN_VEL_Y) > filterParameters_.maxLinVel(1)
    || lastState(State3d::LIN_VEL_Y) < filterParameters_.minLinVel(1)
    || lastState(State3d::LIN_VEL_Z) > filterParameters_.maxLinVel(2)
    || lastState(State3d::LIN_VEL_Z) < filterParameters_.minLinVel(2)) {
        stateSanityMask_ |= 0b100;
    }

    // Check 3D angular velocity
    if (lastState(State3d::ANG_VEL_X) > filterParameters_.maxAngVel(0)
    || lastState(State3d::ANG_VEL_X) < filterParameters_.minAngVel(0)
    || lastState(State3d::ANG_VEL_Y) > filterParameters_.maxAngVel(1)
    || lastState(State3d::ANG_VEL_Y) < filterParameters_.minAngVel(1)
    || lastState(State3d::ANG_VEL_Z) > filterParameters_.maxAngVel(2)
    || lastState(State3d::ANG_VEL_Z) < filterParameters_.minAngVel(2)) {
        stateSanityMask_ |= 0b1000;
    }

    // Check 3D linear acceleration
    if (lastState(State3d::LIN_ACC_X) > filterParameters_.maxLinAcc(0)
    || lastState(State3d::LIN_ACC_X) < filterParameters_.minLinAcc(0)
    || lastState(State3d::LIN_ACC_Y) > filterParameters_.maxLinAcc(1)
    || lastState(State3d::LIN_ACC_Y) < filterParameters_.minLinAcc(1)
    || lastState(State3d::LIN_ACC_Z) > filterParameters_.maxLinAcc(2)
    || lastState(State3d::LIN_ACC_Z) < filterParameters_.minLinAcc(2)) {
        stateSanityMask_ |= 0b10000;
    }

    // Check accelerometer bias
    if (lastState(State3d::ACC_BIAS_X) > filterParameters_.maxAccBias(0)
    || lastState(State3d::ACC_BIAS_X) < filterParameters_.minAccBias(0)
    || lastState(State3d::ACC_BIAS_Y) > filterParameters_.maxAccBias(1)
    || lastState(State3d::ACC_BIAS_Y) < filterParameters_.minAccBias(1)
    || lastState(State3d::ACC_BIAS_Z) > filterParameters_.maxAccBias(2)
    || lastState(State3d::ACC_BIAS_Z) < filterParameters_.minAccBias(2)) {
        stateSanityMask_ |= 0b100000;
    }

    // Check gyroscope bias
    if (lastState(State3d::GYR_BIAS_X) > filterParameters_.maxGyrBias(0)
    || lastState(State3d::GYR_BIAS_X) < filterParameters_.minGyrBias(0)
    || lastState(State3d::GYR_BIAS_Y) > filterParameters_.maxGyrBias(1)
    || lastState(State3d::GYR_BIAS_Y) < filterParameters_.minGyrBias(1)
    || lastState(State3d::GYR_BIAS_Z) > filterParameters_.maxGyrBias(2)
    || lastState(State3d::GYR_BIAS_Z) < filterParameters_.minGyrBias(2)) {
        stateSanityMask_ |= 0b1000000;
    }

    // Check filter covariance matrix
    covSanityMask_ = 0;  // Reset mask
    MatrixX lastCov = states_.back().P();
    if (lastCov(State3d::POS_ERR_X, State3d::POS_ERR_X) > std::pow(filterParameters_.maxPositionStd(0), 2.0)
    || lastCov(State3d::POS_ERR_Y, State3d::POS_ERR_Y) > std::pow(filterParameters_.maxPositionStd(1), 2.0)
    || lastCov(State3d::POS_ERR_Z, State3d::POS_ERR_Z) > std::pow(filterParameters_.maxPositionStd(2), 2.0)) {
        covSanityMask_ |= 0b1;
    }
    if (lastCov(State3d::ATT_ERR_X, State3d::ATT_ERR_X) > std::pow(filterParameters_.maxAttitudeStd(0), 2.0)
    || lastCov(State3d::ATT_ERR_Y, State3d::ATT_ERR_Y) > std::pow(filterParameters_.maxAttitudeStd(1), 2.0)
    || lastCov(State3d::ATT_ERR_Z, State3d::ATT_ERR_Z) > std::pow(filterParameters_.maxAttitudeStd(2), 2.0)) {
        covSanityMask_ |= 0b10;
    }
    if (lastCov(State3d::LIN_VEL_ERR_X, State3d::LIN_VEL_ERR_X) > std::pow(filterParameters_.maxLinVelStd(0), 2.0)
    || lastCov(State3d::LIN_VEL_ERR_Y, State3d::LIN_VEL_ERR_Y) > std::pow(filterParameters_.maxLinVelStd(1), 2.0)
    || lastCov(State3d::LIN_VEL_ERR_Z, State3d::LIN_VEL_ERR_Z) > std::pow(filterParameters_.maxLinVelStd(2), 2.0)) {
        covSanityMask_ |= 0b100;
    }
    if (lastCov(State3d::ANG_VEL_ERR_X, State3d::ANG_VEL_ERR_X) > std::pow(filterParameters_.maxAngVelStd(0), 2.0)
    || lastCov(State3d::ANG_VEL_ERR_Y, State3d::ANG_VEL_ERR_Y) > std::pow(filterParameters_.maxAngVelStd(1), 2.0)
    || lastCov(State3d::ANG_VEL_ERR_Z, State3d::ANG_VEL_ERR_Z) > std::pow(filterParameters_.maxAngVelStd(2), 2.0)) {
        covSanityMask_ |= 0b1000;
    }
    if (lastCov(State3d::LIN_ACC_ERR_X, State3d::LIN_ACC_ERR_X) > std::pow(filterParameters_.maxLinAccStd(0), 2.0)
    || lastCov(State3d::LIN_ACC_ERR_Y, State3d::LIN_ACC_ERR_Y) > std::pow(filterParameters_.maxLinAccStd(1), 2.0)
    || lastCov(State3d::LIN_ACC_ERR_Z, State3d::LIN_ACC_ERR_Z) > std::pow(filterParameters_.maxLinAccStd(2), 2.0)) {
        covSanityMask_ |= 0b10000;
    }
    if (lastCov(State3d::ACC_BIAS_ERR_X, State3d::ACC_BIAS_ERR_X) > std::pow(filterParameters_.maxAccBiasStd(0), 2.0)
    || lastCov(State3d::ACC_BIAS_ERR_Y, State3d::ACC_BIAS_ERR_Y) > std::pow(filterParameters_.maxAccBiasStd(1), 2.0)
    || lastCov(State3d::ACC_BIAS_ERR_Z, State3d::ACC_BIAS_ERR_Z) > std::pow(filterParameters_.maxAccBiasStd(2), 2.0)) {
        covSanityMask_ |= 0b100000;
    }
    if (lastCov(State3d::GYR_BIAS_ERR_X, State3d::GYR_BIAS_ERR_X) > std::pow(filterParameters_.maxGyrBiasStd(0), 2.0)
    || lastCov(State3d::GYR_BIAS_ERR_Y, State3d::GYR_BIAS_ERR_Y) > std::pow(filterParameters_.maxGyrBiasStd(1), 2.0)
    || lastCov(State3d::GYR_BIAS_ERR_Z, State3d::GYR_BIAS_ERR_Z) > std::pow(filterParameters_.maxGyrBiasStd(2), 2.0)) {
        covSanityMask_ |= 0b1000000;
    }

    isSane_ = (stateSanityMask_ > 0) || (covSanityMask_ > 0) ? false : true;

    return isSane_;
}  // checkSanity

void Ekf3d::quaternionToRollPitchYaw(const Quaternion& _q, double* _roll, double* _pitch, double* _yaw) {
  double num = 2.0 * (_q.y() * _q.z() + _q.w() * _q.x());
  double den = 1.0 - 2.0 * (_q.x() * _q.x() + _q.y() * _q.y());
  *_roll = std::atan2(num, den);
  *_pitch = std::asin(2.0 * (_q.w() * _q.y() - _q.x() * _q.z()));
  num = 2.0 * (_q.x() * _q.y() + _q.w() * _q.z());
  den = 1.0 - 2.0 * (_q.y() * _q.y() + _q.z() * _q.z());
  *_yaw = std::atan2(num, den);
}  // quaternionToRollPitchYaw

}  // namespace ekf
