///
/// \file Measurement.hpp
/// \brief Base class for an EKF measurement
/// \author Lionel GENEVE
/// \date 17/11/2018
/// \version 1.0
///

#pragma once

// ---------- Headers ----------
// STD/C++ headers
#include <iostream>
#include <string>
#include <memory>

// EIGEN headers
#include <Eigen/Dense>

// Project headers
#include "FilterUtils.hpp"

///
/// \namespace ekf
/// \brief ekf namespace
///
namespace ekf {

///
/// \class Measurement
/// \brief A base class for an EKF measurement
///
class Measurement {
 public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    typedef std::unique_ptr<Measurement> Ptr;
    typedef std::unique_ptr<const Measurement> ConstPtr;

    /// @brief Construct a new Measurement object
    Measurement() {
//        PRINT_DEBUG("Measurement::Measurement")
        z_ = VectorX::Zero(1);
        R_ = MatrixX::Identity(1,1);
        id_ = -1;
        sensorName_ = "";
        timestamp_ = -1.;
        mahaThresh_ = 1e6;
    }  // Default constructor

    Measurement(const VectorX& _z, const MatrixX& _R,
        const int& _id, const std::string& _name,
        const double& _tStamp, const double& _mahaTh) :
        z_(_z), R_(_R), id_(_id), sensorName_(_name),
        timestamp_(_tStamp), mahaThresh_(_mahaTh) {
        //PRINT_DEBUG("Measurement::Measurement with id= " << id_ << ", sensorName="
        //        << sensorName_ << ", timestamp=" << timestamp_)
    }

    Measurement(const Measurement& _m) = default; ///< Default copy constructor

    virtual ~Measurement() {
//        PRINT_DEBUG("Measurement::~Measurement")
    }  ///< Virtual destructor needed for polymorphism

//    virtual Measurement& operator=(const Measurement& _m);
    virtual Measurement& operator=(const Measurement& _m) = default;

    inline bool operator<(const Measurement& _m) const {
        return (timestamp_ < _m.timestamp());
    }
    inline bool operator>(const Measurement& _m) const {
        return (timestamp_ > _m.timestamp());
    }
    inline bool operator<=(const Measurement& _m) const {
        return (timestamp_ <= _m.timestamp());
    }
    inline bool operator>=(const Measurement& _m) const {
        return (timestamp_ >= _m.timestamp());
    }
    inline bool operator==(const Measurement& _m) const {
        return (timestamp_ == _m.timestamp() && sensorName_ == _m.sensorName());
    }
    inline bool operator!=(const Measurement& _m) const {
        return !(*this == _m);
    }

    friend std::ostream& operator<<(std::ostream& _out, const Measurement& _m) {
        _out << "Measurement " << _m.id() << ", name=";
        _out << _m.sensorName() << ", @ t=" << _m.timestamp();
        _out << ": z=" << _m.z().transpose();
        return _out;
    }

    virtual std::string toString(void) const {
        std::string str;
        str = "Measurement " + std::to_string(id_) + ", name=";
        str += sensorName_ + " @ t=" + std::to_string(timestamp_) + ": z=(";
        for (int i=0;i<z_.size();i++)
            str += std::to_string(z_(i)) + ",";
        str += ")";
        return str;
    }

//    virtual Measurement* clone() const = 0;
//    virtual void destroy(void) = 0;

    // Getter/setter

    /// Measurement vector
    VectorX& z(void) { return z_; }
    const VectorX& z(void) const { return z_; }

    /// Measurement covariance matrix
    MatrixX& R(void) { return R_; }
    const MatrixX& R(void) const { return R_; }

    /// Measurement ID
    int& id(void) { return id_; }
    const int& id(void) const { return id_; }

    /// Sensor name
    std::string& sensorName(void) { return sensorName_; }
    const std::string& sensorName(void) const { return sensorName_; }
    
    /// Timestamp of the measurement in [s]
    double& timestamp(void) { return timestamp_; }
    const double& timestamp(void) const { return timestamp_; }

    /// Mahalanobis threshold associated to the measurement for the outlier test
    double& mahaThresh(void) { return mahaThresh_; }
    const double& mahaThresh(void) const { return mahaThresh_; }

//    virtual VectorX predict(const VectorX& _x) = 0; // Predict the measurement from the current state
//    virtual MatrixX computeJacobian(const VectorX& _state) = 0; // Compute the measurement Jacobian matrix
 private:
 protected:
    VectorX z_;  ///< Measurement vector
    MatrixX R_;  ///< Measurement covariance matrix
    int id_;  ///< Measurement unique ID
    std::string sensorName_;  ///< Name of the corresponding sensor
    double timestamp_;  ///< Time of the measurement [s]
    double mahaThresh_;  ///< Mahalanobis threshold associated to the measurement
};  // class Measurement

}  // namespace ekf
