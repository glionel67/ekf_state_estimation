///
/// \file Ekf.cpp
/// \brief Extended Kalman filter (EKF) base class for state estimation.
/// \author Lionel GENEVE
/// \date 17/11/2018
/// \version 1.0
///

// ---------- Headers ----------
// STD/C++
#include <iostream>
#include <cmath>
#include <cassert>

// Project
#include "Ekf.hpp"

namespace ekf {

template <class T>
Ekf<T>::Ekf(const FilterParameters& _params) {
    //PRINT_DEBUG("Ekf::Ekf")

    // Get parameters
    filterParameters_ = _params;

    PRINT_DEBUG("Ekf::Ekf: nb iters=" << filterParameters_.nIters)
    //assert(filterParameters_.nIters > 0);
    if (filterParameters_.nIters == 0) {
        throw std::invalid_argument("Ekf::init: nIters error: _nIters == 0");
    }

    PRINT_DEBUG("Ekf::Ekf: temporal window=" << filterParameters_.temporalWindow << " s")
    //assert(filterParameters_.temporalWindow > 0 && filterParameters_.temporalWindow < 120);
    if (filterParameters_.temporalWindow <= 0. || filterParameters_.temporalWindow >= 120.) {
        throw std::invalid_argument(
            "Ekf::init: temporal window error: temporalWindow <= 0. || temporalWindow >= 120.");
    }

    size_t n = static_cast<size_t>(std::ceil(filterParameters_.temporalWindow * 100.));
    states_.reserve(n);
    PRINT_DEBUG("Ekf::Ekf: state buffer reserved size=" << n)

    //bufferSize_ = _bufSize;
    //states_.reserve(bufferSize_);

    isInit_ = false;
}

template <class T>
void Ekf<T>::init(const VectorX& _x0, const MatrixX& _P0, double _tStamp) {
    PRINT_DEBUG("Ekf::init")

    assert(!_x0.array().isNaN().any());
    assert(_x0.allFinite());

    assert(!_P0.array().isNaN().any());
    assert(_P0.allFinite());

    PRINT_DEBUG("Ekf::init: creating first state at t=" << _tStamp << " s")
    if (_tStamp < 0.0) {
        throw std::invalid_argument("Ekf::init: timestamp error: _tStamp < 0");
    }

    // Insert first state in the buffer of states
    states_.clear();
    //states_.push_back(T::Ptr(new State(_x0, _P0, _tStamp)));
    states_.push_back(T(_x0, _P0, _tStamp));

    isInit_ = true;
}  // init

template <class T>
bool Ekf<T>::insertMeasurement(const Measurement& _meas, bool _applyCorrection) {
    //PRINT_DEBUG("Ekf::insertMeasurement: " << _meas)

    if (!isInit_) {
        PRINT_DEBUG("Ekf::insertMeasurement: filter not initialized yet")
        return false;
    } else if (states_.size() == 0) {
        PRINT_DEBUG("Ekf::insertMeasurement: no state in the state buffer")
        return false;
    }

    // Assert good measurement
    assert(!_meas.z().array().isNaN().any());
    assert(!_meas.R().array().isNaN().any());
    assert(_meas.z().allFinite());
    assert(_meas.R().allFinite());

    // Find the state with a timestamp just before the measurement
    double tStamp = _meas.timestamp();  // Measurement timestamp [s]
    int idx = static_cast<int>(states_.size()) - 1;  // Index of the last state in the buffer
    //while ((idx > 0) && (states_[idx]->timestamp() > tStamp))
    while ((idx > 0) && (states_[idx].timestamp() > tStamp))
        idx--;

    if (idx < 0) {
        PRINT_DEBUG("Ekf::insertMeasurement: measurement older than the oldest state")
        return false;
    }

    assert(idx >= 0 && idx < states_.size());
    
    // Time difference between previous state and measurement
    //double dt = tStamp - states_[idx]->timestamp();
    double dt = tStamp - states_[idx].timestamp();
    //PRINT_DEBUG("Ekf::insertMeasurement: dt=" << dt)

    if (dt > deltaTimeThresh_) { // No state corresponding to the state[idx] timestamp
        // Measurement corresponding to an old state
        if (idx < states_.size() - 1) {
            // Delta time between measurement and next state
            //double dt2 = states_[idx+1]->timestamp() - tStamp;
            double dt2 = states_[idx + 1].timestamp() - tStamp;
            //PRINT_DEBUG("Ekf::insertMeasurement: dt2=" << dt2)
            // Is the measurement associated to the state at idx+1 ?
            if (dt2 < deltaTimeThresh_) {
                //PRINT_DEBUG("Ekf::insertMeasurement: measurement added to next state")
                //states_[idx+1]->addMeasurement(_meas);
                states_[idx + 1].addMeasurement(_meas);
            }
            else { // No existing state corresponding to the measurement, create a new one
                // Create a new state using the previous state
                //VectorX x = states_[idx]->x();
                //MatrixX P = states_[idx]->P();
                //states_.insert(states_.begin()+idx+1, State::Ptr(new State(x, P, tStamp)));
                //states_[idx+1]->addMeasurement(_meas);
                VectorX x = states_[idx].x();
                MatrixX P = states_[idx].P();
                states_.insert(states_.begin() + idx + 1, T(x, P, tStamp));
                states_[idx + 1].addMeasurement(_meas);
                //PRINT_DEBUG("Ekf::insertMeasurement: inserting a new state")
            }
        }
        else { // Measurement corresponding to a state not yet predicted
            // Create a new state using the latest estimated state
            //VectorX x = states_[idx]->x();
            //MatrixX P = states_[idx]->P();
            VectorX x = states_[idx].x();
            MatrixX P = states_[idx].P();
            // Append the new state in the buffer
            //states_.push_back(State::Ptr(new State(x, P, tStamp)));
            //states_.back()->addMeasurement(_meas);
            states_.push_back(T(x, P, tStamp));
            states_.back().addMeasurement(_meas);
            //PRINT_DEBUG("Ekf::insertMeasurement: appending a new state")
        }
    }
    else if (dt >= 0. && dt <= deltaTimeThresh_) { // The new measurement correspond to an existing state
        //PRINT_DEBUG("Ekf::insertMeasurement: measurement added to previous state")
        //states_[idx]->addMeasurement(_meas);
            states_[idx].addMeasurement(_meas);
    }
    else { // Measurement too old: discard it
        //PRINT_DEBUG("Ekf::insertMeasurement: dt=" << dt << " < 0, measurement too old!")
        return false;
    }

    if (_applyCorrection) {
        // Re-process the states of the buffer from idx
        if (idx > 0)
            idx--;
        
        iterateFromIdx(idx);
    }

    return true;
}  // insertMeasurement

template <class T>
bool Ekf<T>::insertState(const T& _state) {
    //PRINT_DEBUG("Ekf::insertState");
    // Check if the given state is not too close to the last state in the buffer
    double newStateTimestamp = _state.timestamp();  // [s]
    double lastStateTimestamp = states_.back().timestamp();  // [s]
    double deltaTimestamp = newStateTimestamp - lastStateTimestamp;  // [s]
    if (deltaTimestamp > deltaTimeThresh_) {
        // Add new state to the end of the buffer
        states_.push_back(_state);
        // Predict new state
        const size_t nStates = states_.size();
        predict(&(states_[nStates - 2]), &(states_[nStates - 1]));
        return true;
    } else {
        return false;
    }
}  // insertState

template <class T>
void Ekf<T>::iterateFromIdx(const size_t& _idx) {
    recomputeFromIdx(_idx);

    // Check if the buffer is full and delete the oldest states if necessary
    removeOldestState();
}  // iterateFromIdx

template <class T>
void Ekf<T>::removeOldestState(void) {
    //PRINT_DEBUG("Ekf::removeOldestState")
    // Check if the buffer is full and delete the oldest states if necessary
    //double dt = states_.back()->timestamp() - states_.front()->timestamp();
    double dt = states_.back().timestamp() - states_.front().timestamp();
    while (dt > filterParameters_.temporalWindow) {
        //PRINT_DEBUG("Ekf::removeOldestState: last tStamp="
                //<< states_.back()->timestamp() << ", dt=" << dt
        //        << states_.back().timestamp() << ", dt=" << dt
                //<< ", deleting old state @ " << states_.front()->timestamp())
        //        << ", deleting old state @ " << states_.front().timestamp())
        states_.erase(states_.begin());
        //dt = states_.back()->timestamp() - states_.front()->timestamp();
        dt = states_.back().timestamp() - states_.front().timestamp();
    }
}  // removeOldestState

template <class T>
bool Ekf<T>::correct(T* _state) {
    //PRINT_DEBUG("Ekf::correct: for " << *_state)

    // Check if there are measurements associated to the given state
    size_t nMeas = _state->nMeasurements();
    if (nMeas == 0) {
        return false;
    }

    bool ok = true;
    bool okTmp = true;
    for (size_t i = 0; i < nMeas; i++) {  // Loop over all associated measurements
        // Get measurement
        //Measurement* meas = _state->measurements()[i].get();
        Measurement* meas = &(_state->measurements()[i]);

        //okTmp = correctMeasurement(_state, meas);

        int id = meas->id();  // Get measurement ID
        switch (id) {
            case Z_ACCELERO :  // Accelerometer measurement
                okTmp = correctAccelerometer(_state, meas);
                break;
            case Z_GYRO :  // Gyroscope measurement
                okTmp = correctGyroscope(_state, meas);
                break;
            case Z_IMU :  // IMU (accelero + gyro) measurement
                okTmp = correctAccGyr(_state, meas);
                break;
            case Z_MAG :  // Magnetometer measurement
                throw std::runtime_error("Ekf::correct: mag measurement not implemented yet!");
                break;
            case Z_POS :  // Position measurement
                okTmp = correctPosition(_state, meas);
                break;
            case Z_ATT :  // Attitude measurement
                okTmp = correctAttitude(_state, meas);
                break;
            case Z_LIN_VEL :  // Linear velocity measurement
                okTmp = correctLinVel(_state, meas);
                break;
            case Z_ANG_VEL :  // Angular velocity measurement
                okTmp = correctAngVel(_state, meas);
                break;
            case Z_LIN_ANG_VEL : // Linear & angular velocities measurement
                okTmp = correctLinAngVel(_state, meas);
                break;
            case Z_LIN_ACC :  // Acceleration measurement
                okTmp = correctLinAcc(_state, meas);
                break;
            default:
                throw std::runtime_error("Ekf::correct: measurement id not recognized or not implemented!");
                break;
        }  // end case

        ok = ok && okTmp;
    }  // end for

    return ok;
}  // correct

template <class T>
bool Ekf<T>::correctMeasurement(T* /*_state*/, const Measurement* /*_meas*/) {
    /*
    VectorX x = _state->x();  // State vector
    VectorX deltaX = VectorX::Zero(N_ERROR_STATES);  // Error-state vector
    
    // Measurement function Jacobian matrix
    MatrixX H = _meas->computeJacobian(x);

    // Measurement noise (co)variance
    MatrixX R = _meas->R();

    // Measurement prediction from state vector
    VectorX zPred = _meas->predict(x);

    // Innovation vector
    VectorX nu;
    MatrixX PHt;

    // Innovation matrix
    MatrixX S;
    MatrixX Sinv;

    // Kalman gain
    MatrixX K;

    for (size_t n = 0; n < filterParameters_.nIters; n++) {
        // Predict measurement
        zPred = _meas->predict(x);
        assert(zPred.allFinite());
        assert(!zPred.array().isNaN().any());

        // Innovation vector
        nu = _meas->z() - zPred;

        // Innovation matrix
        H = _meas->computeJacobian(x);
        PHt = _state->P() * H.transpose();
        S = H * PHt + R;
        if (!S.allFinite()) {
            std::cerr << "Ekf::correctMeasurement: S=\n" << S << " problem\n";
            return false;
        }

        // Mahalanobis outlier test
        Sinv = S.inverse();
        double maha = nu.transpose() * Sinv * nu;
        //assert(maha >= 0.);
        if (maha < 0.) {
            std::cerr << "Ekf::correctMeasurement: maha=" << maha << " < 0\n";
            return false;
        }
        //maha = std::sqrt(maha);
        if (maha > _meas->mahaThresh()) {
            std::cerr << "Ekf::correctMeasurement: maha=" << maha << " > thresh=" << _meas->mahaThresh() << std::endl;
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
        deltaX = K * nu;  // Error-state vector
        assert(deltaX.allFinite());
        assert(!deltaX.array().isNaN().any());
        addErrorState(x, deltaX); // x += deltaX;

        // Stop if no further improvement
        if (deltaX.norm() < 1e-6)
            break;
    }  // Iterated correction

    // State vector correction
    assert(x.allFinite());
    assert(!x.array().isNaN().any());
    _state->x() = x;
    
    // Covariance matrix correction
    //MatrixX I_ = MatrixX::Identity(N_ERROR_STATES, N_ERROR_STATES);
    MatrixX IKH = I_ - K * H;
    _state->P() = IKH * _state->P() * IKH.transpose() + K * R * K.transpose();
    errorReset(_state->P(), deltaX);
    
    // Check state vector and covariance matrix
    _state->checkStateAndCovariance();
    */
    
    return true;
}  // correctMeasurement

template <class T>
void Ekf<T>::recomputeFromIdx(const size_t& _idx) {
    //PRINT_DEBUG("Ekf::recomputeFromIdx: idx=" << _idx)

    //assert(_idx >= 0 && _idx < states_.size() - 1);
    if (_idx >= states_.size() - 1) {
        throw std::invalid_argument("Ekf::recomputeFromIdx: index error: _idx >= states_.size()-1");
    }

    size_t currIdx = 0, nextIdx = 0;
    for (size_t i = _idx; i < states_.size() - 1; i++) {
        currIdx = i;
        nextIdx = i + 1;
        // Prediction from current (idx) to next (idx+1) state
        //PRINT_DEBUG("Ekf::recomputeFromIdx: prediction from " << currIdx << " to " << nextIdx)
        //predict(states_[currIdx].get(), states_[nextIdx].get());
        predict(&(states_[currIdx]), &(states_[nextIdx]));

        // Check if there are measurements associated to the next state (idx+1)
        //if (states_[nextIdx]->measurements().size() > 0) {
        if (states_[nextIdx].measurements().size() > 0) {
            //PRINT_DEBUG("Ekf::recomputeFromIdx: correction for " << nextIdx)
            //correct(states_[nextIdx].get());
            correct(&(states_[nextIdx]));
        }
    }
}  // recomputeFromIdx

template <class T>
T* Ekf<T>::getLastState(void) {
    if (!states_.empty()) {
        //return states_.back().get();
        return &(states_.back());
    }
    else {
        return nullptr;
    }
}  // getLastState

template <class T>
void Ekf<T>::reset(void) {
    //PRINT_DEBUG("Ekf::reset");
    states_.clear();  // Remove/clear all the state
    isInit_ = false;  // Force re-init
    // Reallocate memory for state buffer
    size_t n = static_cast<size_t>(std::ceil(filterParameters_.temporalWindow * 100.));
    states_.reserve(n);
}  // reset

}  // namespace ekf

#include "State2d.hpp"
#include "State3d.hpp"
namespace ekf {
  template class Ekf<State2d>;
  template class Ekf<State3d>;
}
