///
/// \file FilterUtils.hpp
/// \brief Filter utility defines and functions
/// \author Lionel GENEVE
/// \date 17/11/2018
/// \version 1.0
///

#pragma once

// ---------- Headers ----------
// STD/C++ headers
#include <iostream>
#include <cmath>

// EIGEN headers
#include <Eigen/Dense>

// ---------- Defines ----------
#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

#define DEG_2_RAD (M_PI / 180.0)
#define RAD_2_DEG (180.0 / M_PI)

#define DOUBLE_PRECISION  // Double 64 precision
//#define FLOAT_PRECISION  // Uncomment to activate float 32 precision

//#define FILTER_2D  // Uncomment to activate the 2D filter version

#define DEBUG  // Uncomment to disable debug mode

#ifdef DEBUG
    #define PRINT_DEBUG(msg) std::cerr << msg << std::endl;
#else
    #define PRINT_DEBUG(msg)
#endif

///
/// \namespace ekf
/// \brief ekf namespace
///
namespace ekf {

#ifdef FLOAT_PRECISION
typedef Eigen::VectorXf VectorX;
typedef Eigen::Vector2f Vector2;
typedef Eigen::Vector3f Vector3;
typedef Eigen::Vector4f Vector4;

typedef Eigen::MatrixXf MatrixX;
typedef Eigen::Matrix22f Matrix22;
typedef Eigen::Matrix33f Matrix33;
typedef Eigen::Matrix44f Matrix44;

typedef Eigen::Quaternionf Quaternion;
#endif

#ifdef DOUBLE_PRECISION
typedef Eigen::VectorXd VectorX;
typedef Eigen::Vector2d Vector2;
typedef Eigen::Vector3d Vector3;
typedef Eigen::Vector4d Vector4;

typedef Eigen::MatrixXd MatrixX;
typedef Eigen::Matrix2d Matrix2;
typedef Eigen::Matrix3d Matrix3;
typedef Eigen::Matrix4d Matrix4;

typedef Eigen::Quaterniond Quaternion;
#endif

///
/// \enum MeasurementIdx_t
/// \brief Index (IDs) of the different measurements fused in the filter
///
typedef enum {
    Z_ACCELERO = 0,
    Z_GYRO,  // 1
    Z_IMU,  // 2
    Z_MAG,  // 3
    Z_POS,  // 4
    Z_ATT,  // 5
    Z_LIN_VEL,  // 5
    Z_ANG_VEL,  // 6
    Z_LIN_ANG_VEL,  // 7
    Z_LIN_ACC,  // 8
    Z_YAW,  // 9
    Z_GPS_POS,  // 10
    Z_GPS_YAW,  // 11
    Z_GPS_VEL,  // 12
    N_MEASUREMENTS  // 13
} MeasurementIdx_e;

///
/// \fn constrainValue
/// \brief Constrain a value between a min and a max value and check for NaN and Inf
/// \param _val: (in/out) value to check
/// \param _min: minimum value
/// \param _max: maximal value
///
inline void constrainValue(double& _val, const double& _min, const double& _max);

/**
 * @brief Wrap angle in (-pi, pi]
 * @fn piToPi
 * @param _angle angle to wrap in [rad]
 * @return double wrapped angle in [rad]
 */
double piToPi(double _angle);

} // namespace ekf
