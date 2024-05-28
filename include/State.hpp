///
/// \file State.hpp
/// \brief Base class for an EKF state
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
#include "FilterUtils.hpp"
#include "Measurement.hpp"

// ---------- Defines ----------

///
/// \namespace ekf
/// \brief ekf namespace
///
namespace ekf {

///
/// \class State
/// \brief A base class for an EKF state
///
class State {
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  typedef std::unique_ptr<State> Ptr;
  typedef std::unique_ptr<const State> ConstPtr;

  State() {
    //PRINT_DEBUG("State::State")
    timestamp_ = -1.;
    measurements_.clear();
    measurements_.reserve(nMeasAllocSize);
  }  ///< Default constructor

  State(double _tStamp) : timestamp_(_tStamp) {
    //PRINT_DEBUG("State::State @ t=" << timestamp_)
    measurements_.clear();
    measurements_.reserve(nMeasAllocSize);
  }

  State(const VectorX& _x, const MatrixX& _P, const double& _tStamp) : 
            x_(_x), P_(_P), timestamp_(_tStamp) {
    //PRINT_DEBUG("State::State @ t=" << timestamp_ << ", with x=" << x_.transpose())
    measurements_.clear();
    measurements_.reserve(nMeasAllocSize);
  }

  State(const State& _s) = default;  ///< Default copy constructor

  virtual ~State() {
//    PRINT_DEBUG("State::~State")
  }  ///< Default destructor

  State& operator=(const State& _s) = default;

    inline bool operator<(const State& _s) const {
        return (timestamp_ < _s.timestamp());
    }
    inline bool operator>(const State& _s) const {
        return (timestamp_ > _s.timestamp());
    }
    inline bool operator<=(const State& _s) const {
        return (timestamp_ <= _s.timestamp());
    }
    inline bool operator>=(const State& _s) const {
        return (timestamp_ >= _s.timestamp());
    }
    inline bool operator==(const State& _s) const {
        return (timestamp_ == _s.timestamp());
    }
    inline bool operator!=(const State& _s) const {
        return !(*this == _s);
    }

    friend std::ostream& operator<<(std::ostream& _out, const State& _s) {
        _out << "State @ t=" << _s.timestamp();
        _out << " with nMeas=" << _s.measurements().size();
        _out << ", x=" << _s.x().transpose();
        //_out << ", P=\n" << _s.P();
        return _out;
    }

    std::string toString(void) const {
        std::string str;
        str = "State @ t=" + std::to_string(timestamp_) + " with nMeas=";
        str += std::to_string(measurements_.size()) + ", x=(";
        for (int i=0;i<x_.size();i++)
            str += std::to_string(x_(i)) + ",";
        str += ")";
        return str;
    }

    ///
    /// \fn addErrorState
    /// \brief Add the error state to the state vector
    /// \param _dx: Error-state vector to inject into the state vector
    /// 
    void addErrorState(const VectorX& _dx) {
        // TODO
    }

    ///
    /// \fn addMeasurement
    /// \brief Associated a given measurement to the state
    /// 
    void addMeasurement(const Measurement& _m) {
        //measurements_.push_back(Measurement::Ptr(new Measurement(_m)));
        measurements_.push_back(_m);
    }

    //void addErrors(const VectorX& _dx);
    //void resetCovariance(const VectorX& _dx);
    //void reset(void);

  ///
  /// \fn checkStateVector
  /// \brief Check the state vector from inf and nan values
  /// 
  virtual void checkStateVector(void) {
    //virtual void checkStateVector(const VectorX& xMin, const VectorX& xMax) {
        //assert(!x_.array().isNaN().any());
        //assert(x_.allFinite());
        for (int i = 0; i < x_.size(); i++) {
            if (std::isnan(x_(i))) {
                std::cerr << *this << " --> isnan at index [" << i << "]!\n";
                x_(i) = 0.;
            }
            if (std::isinf(x_(i))) {
                std::cerr << *this << " --> isinf at index [" << i << "]!\n";
                x_(i) = 0.;
            }
            //constrainValue(x_(i), xMin(i), xMax(i));
        }
  }  // checkStateVector
    
  ///
  /// \fn checkCovarianceMatrix
  /// \brief Check the covariance matrix from inf and nan values
  ///
  virtual void checkCovarianceMatrix(void) {
    //virtual void checkCovarianceMatrix(double Pmin, double Pmax) {
        //assert(!P_.array().isNaN().any());
        //assert(P_.allFinite());
        for (int i = 0; i < P_.rows(); i++) {
            for (int j = 0; j < P_.cols(); j++) {
                if (std::isnan(P_(i, j))) {
                    std::cerr << *this << " --> isnan at index [" << i << ", " << j << "]!\n";
                    P_(i,j) = 1.;
                }
                if (std::isinf(P_(i, j))) {
                    std::cerr << *this << " --> isinf at index [" << i << ", " << j << "]!\n";
                    P_(i,j) = 1.;
                }
                //constrainValue(P_(i,j), Pmin, Pmax);
            }
        }
  }  // checkCovarianceMatrix

  ///
  /// \fn checkStateAndCovariance
  /// \brief Check the state vector and covariance matrix from inf and nan values
  ///
  virtual void checkStateAndCovariance(void) {
  //virtual void checkStateAndCovariance(const VectorX& xMin, const VectorX& xMax, double Pmin, double Pmax) {
    //assert(!x_.array().isNaN().any());
    //assert(x_.allFinite());
    //assert(!P_.array().isNaN().any());
    //assert(P_.allFinite());

    for (int i = 0; i < P_.rows(); i++) {
      if (std::isnan(x_(i)) || std::isinf(x_(i))) {
        std::cerr << *this << " --> isnan/isinf at index [" << i << "]!\n";
        x_(i) = 0.;
      }
      //constrainValue(x_(i), xMin(i), xMax(i));
      for (int j = 0; j < P_.cols(); j++) {
        if (std::isnan(P_(i, j)) || std::isinf(P_(i, j))) {
          std::cerr << *this << " --> isnan/isinf at index [" << i << ", " << j << "]!\n";
          P_(i, j) = 1.;
        }
        //constrainValue(P_(i,j), Pmin, Pmax);
      }
    }
  }  // checkStateAndCovariance

  ///
  /// \fn enforceCovMatSymmetry
  /// \brief Enforce the covariance matrix to be symmetric
  /// 
  void enforceCovMatSymmetry(void) {
    /*double p = 0.0;

    for (int i = 0; i < P_.rows(); i++) {
        for (int j = 0; j < P_.cols(); j++) {
            p = .5 * (P_(i, j) + P_(j, i));
            P_(i, j) = P_(j, i) = p;
        }
    }*/

    P_ = 0.5 * (P_ + P_.transpose());
  }  // enforceCovMatSymmetry

  size_t nMeasurements(void) const { return measurements_.size(); }  // Return number of measurements associated to the state

  // Getter/setter
  VectorX& x(void) { return x_; }
  const VectorX& x(void) const { return x_; }

  MatrixX& P(void) { return P_; }
  const MatrixX& P(void) const { return P_; }
    
  double& timestamp(void) { return timestamp_; }
  const double& timestamp(void) const { return timestamp_; }

//  std::vector<Measurement::Ptr>& measurements(void) { return measurements_; }
//  const std::vector<Measurement::Ptr>& measurements(void) const { return measurements_; }

  std::vector<Measurement>& measurements(void) { return measurements_; }
  const std::vector<Measurement>& measurements(void) const { return measurements_; }
 private:
 protected:
  VectorX x_;  ///< State vector
  MatrixX P_;  ///< Covariance matrix associated to the state vector
  double timestamp_;  ///< Timestamp of the state [s]
//  std::vector<Measurement::Ptr> measurements_;  ///< Vector of measurements associated to the state
  std::vector<Measurement> measurements_;  ///< Vector of measurements associated to the state
//  int nbStates_;  ///< Size of the state vector
//  int nbErrorStates_;  ///< Size of the error-state vector

  static constexpr size_t nMeasAllocSize = 3;
};  // class State

}  // namespace ekf
