///
/// \file Ekf.hpp
/// \brief Extended Kalman filter (EKF) base class for state estimation.
/// \author Lionel GENEVE
/// \date 17/11/2018
/// \version 1.0
///

#pragma once

// ---------- Headers ----------
// STD/C++ headers
#include <vector>
#include <map>

// EIGEN headers
#include <Eigen/Core>
#include <Eigen/Dense>

// Project headers
#include "FilterParameters.hpp"
#include "FilterUtils.hpp"
#include "Measurement.hpp"
//#include "RingBuffer.hpp"
//#include "State.hpp"

///
/// \namespace ekf
/// \brief ekf namespace
///
namespace ekf {

///
/// \class Ekf
/// \brief Implementation of a 3D state estimation filter using an Extended Kalman filter (EKF)
///
template <class T>
class Ekf {
 protected:
  //typedef std::map<double, T> StateBuffer;
  //typedef RingBuffer<T> StateBuffer;
  //typedef std::vector<T::Ptr> StateBuffer;
  typedef std::vector<T> StateBuffer;
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  Ekf() = delete;  ///< Delete default constructor

  ///
  /// \fn Ekf
  /// \brief Ekf base class constructor with parameters
  /// \param _params: filter parameters
  ///
  Ekf(const FilterParameters& _params);

  /**
   * @brief Copy constructor
   * @param _ekf 
   */
  Ekf(const Ekf& _ekf) = default;  ///< Default copy constructor

  /**
   * @brief Destroy the Ekf object
   * 
   */
  virtual ~Ekf() = default;  ///< Default destructor

  /**
   * @brief Default assignement operator
   * @param _ekf 
   * @return Ekf& 
   */
  Ekf& operator=(const Ekf& _ekf) = default;

  ///
  /// \fn reset
  /// \brief Reset the filter
  ///
  void reset(void);

  ///
  /// \fn init
  /// \brief Initialize the EKF with first state
  /// \param _x0: initial state vector
  /// \param _P0: initial covariance matrix
  /// \param _tStamp: timestamp of the initial state [s]
  ///
  virtual void init(const VectorX& _x0, const MatrixX& _P0, double _tStamp);

  ///
  /// \fn insertMeasurement
  /// \brief Insert a new measurement in the filter (and apply correction)
  /// \param _meas: measurement to insert in the filter
  /// \param _applyCorrection: if the correction step must be performed
  /// \return true if new measurement was fused in the filter, false else
  ///
  bool insertMeasurement(const Measurement& _meas, bool _applyCorrection = true);

  ///
  /// \fn insertState
  /// \brief Insert a new state in the filter (and apply prediction)
  /// \param _state: state to insert in the filter
  /// \return true if new state was added to the state buffer, false else
  ///
  bool insertState(const T& _state);

  ///
  /// \fn iterateFromIdx
  /// \brief Apply the prediction and correction steps from given index to last (newest) state
  /// \param _idx: index of the state buffer from which to start recompute the filter
  ///
  void iterateFromIdx(const size_t& _idx);

  ///
  /// \fn getLastState
  /// \brief Return a pointer to the last (newest) estimated state
  ///
  T* getLastState(void);

  ///
  /// \fn checkSanity
  /// \brief Check filter sanity
  ///
  virtual bool checkSanity(void) {return 0;}

  // Getters/setters
  StateBuffer& states(void) { return states_; }
  const StateBuffer& states(void) const { return states_; }

  FilterParameters& filterParameters(void) { return filterParameters_; }
  const FilterParameters& filterParameters(void) const { return filterParameters_; }

  bool& isInit(void) { return isInit_; }
  const bool& isInit(void) const { return isInit_; }

  uint8_t& stateSanityMask(void) { return stateSanityMask_; }
  const uint8_t& stateSanityMask(void) const { return stateSanityMask_; }

  uint8_t& covSanityMask(void) { return covSanityMask_; }
  const uint8_t& covSanityMask(void) const { return covSanityMask_; }

  ///
  /// \fn recomputeFromIdx
  /// \brief Recompute the prediction and correction step from the given index to last (newest) state
  /// \param _idx: index of the state buffer from which to start recompute the filter
  ///
  void recomputeFromIdx(const size_t& _idx);

  ///
  /// \fn removeOldestState
  /// \brief Remove states older than the time window from the state buffer
  ///
  void removeOldestState(void);

  ///
  /// \fn predict
  /// \brief Predict the next state from the current state
  /// \param _currState: current state
  /// \param _nextState: next state to predict
  ///
  virtual void predict(const T* _currState, T* _nextState) {}

  ///
  /// \fn correct
  /// \brief Apply the correction step of the EKF on the state
  /// \param _state: state for which to apply the correction
  /// \return bool: true if the correction succeeded, else false
  ///
  bool correct(T* _state);

  ///
  /// \fn correctMeasurement
  /// \brief Apply the correction step on the state with the given measurement
  /// \param _state: pointer to the state for which to apply the correction
  /// \param _meas: pointer to the measurement used to make the correction
  /// \return bool: true if the correction succeeded, else false
  ///
  virtual bool correctMeasurement(T* _state, const Measurement* _meas);
  
  virtual bool correctAccelerometer(T* _state, const Measurement* _meas)  {return 0;}
  virtual bool correctGyroscope(T* _state, const Measurement* _meas)  {return 0;}
  virtual bool correctAccGyr(T* _state, const Measurement* _meas)  {return 0;}
  virtual bool correctLinVel(T* _state, const Measurement* _meas)  {return 0;}
  virtual bool correctAngVel(T* _state, const Measurement* _meas)  {return 0;}
  virtual bool correctLinAngVel(T* _state, const Measurement* _meas)  {return 0;}
  virtual bool correctLinAcc(T* _state, const Measurement* _meas)  {return 0;}
  virtual bool correctPosition(T* _state, const Measurement* _meas)  {return 0;}
  virtual bool correctAttitude(T* _state, const Measurement* _meas)  {return 0;}
  //virtual bool correctGpsPosition(T* _state, const Measurement* _meas)  {return 0;}
  //virtual bool correctGpsSpeed(T* _state, const Measurement* _meas)  {return 0;}

  ///
  /// \fn addErrorState
  /// \brief Add the error-state to the state vector
  /// \param _x: state vector to update
  /// \param _dx: error-state vector to add
  ///
  virtual void addErrorState(VectorX& _x, VectorX& _dx) {}

  ///
  /// \fn errorReset
  /// \brief Reset the error from the covariance matrix
  /// \param _P: covariance matrix to update
  /// \param _dx: error-state vector used for the reset
  ///
  virtual void errorReset(MatrixX& _P, VectorX& _dx) {}

 private:
 protected:
  StateBuffer states_;  ///< Buffer of states
  FilterParameters filterParameters_;  ///< Parameters of the filter
  bool isInit_;  ///< If the filter is initalized with a state vector and covariance matrix
  MatrixX I_;  ///< Identity matrix
  MatrixX F_;  ///< State transition matrix
  bool isSane_{true};  ///< If the filter is running correctly
  uint8_t stateSanityMask_{0};  ///< Bitfield mask to track which state is a problem in the state vector
  uint8_t covSanityMask_{0};  ///< Bitfield mask to track which state is a problem in the covariance matrix

  static constexpr double deltaTimeThresh_ = 5e-3;  ///< Threshold used to associated a measurement with a state [s]
};  // class Ekf

}  // namespace ekf

//#include "Ekf.cpp"