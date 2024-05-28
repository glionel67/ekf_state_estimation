///
/// \file test_ekf.cpp
/// \brief Test function of the EKF
/// \author Lionel GENEVE
/// \date 30/03/2022
/// \version 1.0
///

// ---------- Headers ----------
// STD
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <queue>
#include <utility>
#include <chrono>
#include <bitset>

// EIGEN
#include <Eigen/Dense>

// Project
#include "GpsConversion.hpp"
#include "FilterUtils.hpp"
#include "Measurement.hpp"
#include "State2d.hpp"
#include "State3d.hpp"
#include "Ekf2d.hpp"
#include "Ekf3d.hpp"

using namespace ekf;
//using namespace std::chrono;
namespace chr = std::chrono;

// Useful functions
bool loadDataset(const std::string& _filename, std::vector<std::vector<double>>& _dataset) {
    // Clear dataset
    _dataset.clear();
    _dataset.reserve(10000);

    // Open file
    std::ifstream file(_filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open dataset file: " << _filename << std::endl;
        return false;
    }

    // Parse file line by line
    size_t nLines = 0;
    std::string line;
    while (std::getline(file, line)) {
        std::stringstream lineStream(line);
        std::string cell;
        std::vector<double> result;
        result.reserve(30);

        while (std::getline(lineStream, cell, ';')) {
            result.push_back(std::atof(cell.c_str()));
        }

        _dataset.push_back(result);
        nLines++;
    }

    std::cout << "Parsed dataset with nLines=" << nLines << std::endl;
    
    return true;
}  // loadDataset

void createFilterParameters(FilterParameters* filterParams) {
  filterParams->positionProcessNoiseStd = 1e-4 / 3.;
  filterParams->attitudeProcessNoiseStd = (.1 * M_PI / 180.) / 3.;
  filterParams->linVelProcessNoiseStd = 1e-2;
  filterParams->linAccProcessNoiseStd = 1e-2;
  filterParams->angVelProcessNoiseStd = 1e-3;
  //filterParams->accBiasNoiseStd = 1e-5;  // RT2000
  //filterParams->gyrBiasNoiseStd = 1e-6;  // RT2000
  filterParams->accBiasNoiseStd = 1e-4;  // PX4
  filterParams->gyrBiasNoiseStd = 1e-5;  // PX4

  filterParams->nIters = 1;
  filterParams->temporalWindow = 1.0;
  filterParams->minCovMatrix = 1e-5;
  filterParams->maxCovMatrix = 1e3;
  filterParams->minPosition << -10e3, -10e3, -10e3;
  filterParams->maxPosition << 10e3, 10e3, 10e3;
  filterParams->minAttitude << -M_PI, -M_PI, -M_PI;
  filterParams->maxAttitude << M_PI, M_PI, M_PI;
  filterParams->minLinVel << -30.0, -10.0, -10.0;
  filterParams->maxLinVel << 300., 10.0, 10.0;
  filterParams->minAngVel << -7.0, -7.0, -7.0; 
  filterParams->maxAngVel << 7.0, 7.0, 7.0;
  filterParams->minLinAcc << -10.0, -4.0, -4.0;
  filterParams->maxLinAcc << 8.0, 4.0, 4.0;
  filterParams->minAccBias << -10.0, -4.0, -4.0;
  filterParams->maxAccBias << 8.0, 4.0, 4.0;
  filterParams->minGyrBias << -7.0, -7.0, -7.0;
  filterParams->maxGyrBias << 7.0, 7.0, 7.0;

  filterParams->maxPositionStd << 10.0, 10.0, 10.0;
  filterParams->maxAttitudeStd << 0.5, 0.5, 0.5;
  filterParams->maxLinVelStd << 1.0, 1.0, 1.0;
  filterParams->maxAngVelStd << 0.2, 0.2, 0.2;
  filterParams->maxLinAccStd << 0.5, 0.5, 0.5;
  filterParams->maxAccBiasStd << 0.3, 0.3, 0.3;
  filterParams->maxGyrBiasStd << 0.1, 0.1, 0.1;

  filterParams->gravityVector << 0.0, 0.0, -9.81;
}  // createFilterParameters

std::string row2str(const std::vector<double>& _row) {
    std::string str;
    if (_row.size() <= 0) {
        return str;
    }

    for (const auto& val : _row) {
        str += std::to_string(val) + ", ";  
    }
    //str += "\n";

    return str;
}  // row2str

std::string stateToString(const VectorX& _x, const MatrixX& _P) {
    std::stringstream ss;

    //std::streamsize p = ss.precision();
    ss.precision(12);

    ss << _x(0);
    for (int i = 1; i < _x.size(); i++) {
        ss << ";" << _x(i);
    }

    for (int i = 0; i < _P.rows(); i++) {
        ss << ";" << _P(i, i);
    }

    ss << "\n";

    return ss.str();
}  // stateToString

void test_ekf3(const std::vector<std::vector<double>>& dataset, const std::string& resultFilename) {
    // Output file
    std::ofstream outFile(resultFilename);

    // Create filter parameters
    FilterParameters filterParams;
    createFilterParameters(&filterParams);

    // Measurement noise covariance/std
    // --- RT2000
    const double acceleroNoiseStd = 4.0e-2; // [m/s2]
    const double gyroNoiseStd = 3.0e-3; // [rad/s]
    // --- PX4
    //const double acceleroNoiseStd = 2.0e-2; // [m/s2]
    //const double gyroNoiseStd = 2.0e-4; // [rad/s]

    const double linVelNoiseStdXY = 0.3; // [m/s]
    const double linVelNoiseStdZ = 0.1; // [m/s]
    const double angVelNoiseStdXY = 0.05; // [rad/s]
    const double angVelNoiseStdZ = 0.03; // [rad/s]

    // Create estimation filter
    Ekf3d* ekf = new Ekf3d(filterParams);

    const int N_STATES = State3d::N_STATES;
    const int N_ERROR_STATES = State3d::N_ERROR_STATES;

    // Initialize first state of the estimation filter
    VectorX x0 = VectorX::Zero(N_STATES); // vector length = 22
    x0.segment(State3d::POSITION_X, 3) << 0.0, 0.0, 0.0; // 3D position
    x0.segment(State3d::QUATERNION_W, 4) << 1.0, 0.0, 0.0, 0.0; // 3D attitute (quaternion)
    x0.segment(State3d::LIN_VEL_X, 3) << 0.0, 0.0, 0.0; // 3D linear velocity
    x0.segment(State3d::ANG_VEL_X, 3) << 0.0, 0.0, 0.0; // 3D angular velocity
    x0.segment(State3d::LIN_ACC_X, 3) << 0.0, 0.0, 0.0; // 3D linear acceleration
    x0.segment(State3d::ACC_BIAS_X, 3) << 0.0, 0.0, 0.0; // Accelero bias
    x0.segment(State3d::GYR_BIAS_X, 3) << 0.0, 0.0, 0.0; // Gyro bias

    MatrixX P0 = MatrixX::Identity(N_ERROR_STATES, N_ERROR_STATES);
    P0.block<3, 3>(State3d::POS_ERR_X, State3d::POS_ERR_X) = std::pow(0.3, 2.0) * Matrix3::Identity();
    P0.block<3, 3>(State3d::ATT_ERR_X, State3d::ATT_ERR_X) = std::pow(0.03, 2.0) * Matrix3::Identity();
    P0.block<3, 3>(State3d::LIN_VEL_ERR_X, State3d::LIN_VEL_ERR_X) = std::pow(0.1, 2.0) * Matrix3::Identity();
    P0.block<3, 3>(State3d::ANG_VEL_ERR_X, State3d::ANG_VEL_ERR_X) = std::pow(0.01, 2.0) * Matrix3::Identity();
    P0.block<3, 3>(State3d::LIN_ACC_ERR_X, State3d::LIN_ACC_ERR_X) = std::pow(0.01, 2.0) * Matrix3::Identity();
    P0.block<3, 3>(State3d::ACC_BIAS_ERR_X, State3d::ACC_BIAS_ERR_X) = std::pow(0.15, 2.0) * Matrix3::Identity();
    P0.block<3, 3>(State3d::GYR_BIAS_ERR_X, State3d::GYR_BIAS_ERR_X) = std::pow(0.05, 2.0) * Matrix3::Identity();

    const double t0 = dataset[0][0];  // [s]

    // Initialize filter
    ekf->init(x0, P0, t0);
    outFile << stateToString(ekf->getLastState()->x(), ekf->getLastState()->P());  // Add to result

    // GPS initialization
    double gpsLat0Rad = dataset[0][19] * DEG_2_RAD;
    double gpsLon0Rad = dataset[0][20] * DEG_2_RAD;
    double gpsAlt0M = dataset[0][21];
    double ecefX = 0.0, ecefY = 0.0, ecefZ = 0.0;
    double east = 0.0, north = 0.0, up = 0.0;
    llaToEcef(gpsLat0Rad, gpsLon0Rad, gpsAlt0M, &ecefX, &ecefY, &ecefZ);
    ecefToEnu(ecefX, ecefY, ecefZ, gpsLat0Rad, gpsLon0Rad, gpsAlt0M, &east, &north, &up);

    // Initialize useful variables
    double previousTimestamp = t0;  // [s]
    double currentTimestamp = t0;  // [s]
    const double timestampThresh = 0.01;  // [s]
    bool vehicleAtRest = false;  // If the vehicle is stopped and not moving

    const double wheelRadiusM = 0.3575;  // [m]
    const double dFrontAxleToCog = 1.135;  // [m]
    const double dRearAxleToCog = 1.115;  // [m]
    const double wheelBaseM = dFrontAxleToCog + dRearAxleToCog;  // [m]
    const double rpm2rads = M_PI / 30.0;

    // Timing variables
    chr::high_resolution_clock::time_point tStart, tStop;
    double tSpanSec = 0.0;  // [s]

    // Main temporal loop over the dataset
    //bool isRunning = true;
    const size_t nRows = dataset.size();
    size_t n = 1;
    //while (isRunning) {
    while (n < nRows) {
        currentTimestamp = dataset[n][0];  // [s]
        if (0 == n % 100)
            std::cout << "Iteration " << n << ": time=" << currentTimestamp << " s" << std::endl;

        // Create new state and add it to the filter
        double deltaTimestamp = currentTimestamp - ekf->getLastState()->timestamp();
        //std::cout << "Last state timestamp=" << ekf->getLastState()->timestamp() << " s" << std::endl;
        if (deltaTimestamp >= timestampThresh) {
            //std::cerr << "Creating a new state..." << std::endl;
            State3d state(currentTimestamp);
            //tStart = chr::high_resolution_clock::now();
            ekf->insertState(state);
            //tStop = chr::high_resolution_clock::now();
            //tSpanSec = chr::duration_cast<chr::duration<double>>(tStop - tStart).count();
            //std::cout << "Took " << tSpanSec << " s to insert a new state" << std::endl;
        }

        // Check for new measurements

        // IMU measurement
        double zax = dataset[n][1], zay = dataset[n][2], zaz = dataset[n][3];  // RT2000 accelero
        double zgx = dataset[n][4], zgy = dataset[n][5], zgz = dataset[n][6];  // RT2000 gyro

        double zax2 = dataset[n][7], zay2 = dataset[n][8], zaz2 = dataset[n][9];  // PX4 accelero
        double zgx2 = dataset[n][10], zgy2 = dataset[n][11], zgz2 = dataset[n][12];  // PX4 gyro

        VectorX z_imu = VectorX::Zero(6);
        //z_imu << zax, zay, zaz, zgx, zgy, zgz;
        z_imu << zax2, zay2, zaz2, zgx2, zgy2, zgz2;
        MatrixX R_imu = MatrixX::Identity(6, 6);
        R_imu(0, 0) = R_imu(1, 1) = R_imu(2, 2) = std::pow(acceleroNoiseStd, 2.0);
        R_imu(3, 3) = R_imu(4, 4) = R_imu(5, 5) = std::pow(gyroNoiseStd, 2.0);
        int id_imu = Z_IMU;
        std::string name_imu = "imu";
        double ts_imu = currentTimestamp;
        double mahaThresh_imu = 12.59;  // chi2(0.05, 6)
        Measurement meas_imu(z_imu, R_imu, id_imu, name_imu, ts_imu, mahaThresh_imu);
        ekf->insertMeasurement(meas_imu, false);

        // Vehicle measurements
        int transmissionMode = static_cast<int>(dataset[n][24]);  // 0: sna, 1: park, 2: reverse, 3: neutral, 4: forward
        double vehicleSpeedKmh = dataset[n][25];  // [km/h]

        vehicleAtRest = (((1 == transmissionMode) || (3 == transmissionMode)) && (std::abs(vehicleSpeedKmh) <= 0.1)) ? true : false;

        if (vehicleAtRest) {  // Vehicle at rest: ZUP, ZUV, ZUA
            // Zero Update Velocities (ZUV): force linear and angular velocities to zero
            VectorX z_zuv = VectorX::Zero(6);
            MatrixX R_zuv = MatrixX::Identity(6, 6);
            R_zuv(0, 0) = R_zuv(1, 1) = R_zuv(2, 2) = std::pow(0.01, 2.0);  // Around 0.1 km/h 3*std
            R_zuv(3, 3) = R_zuv(4, 4) = R_zuv(5, 5) = std::pow(0.001, 2.0);  // Around 0.2 deg/s 3*std
            int id_zuv = Z_LIN_ANG_VEL;
            std::string name_zuv = "veh_zuv";
            double ts_zuv = currentTimestamp;
            double mahaThresh_zuv = 12.59;  // chi2(0.05, 6)
            Measurement meas_zuv(z_zuv, R_zuv, id_zuv, name_zuv, ts_zuv, mahaThresh_zuv);
            ekf->insertMeasurement(meas_zuv, false);

            // Zero Update Acceleration (ZUA): force linear acceleration to zero
            VectorX z_zua = VectorX::Zero(3);
            MatrixX R_zua = MatrixX::Identity(3, 3) * std::pow(0.0033, 2.0);  // Around 0.01 m/s2 3*std
            int id_zua = Z_LIN_ACC;
            std::string name_zua = "veh_zua";
            double ts_zua = currentTimestamp;
            double mahaThresh_zua = 7.81;  // chi2(0.05, 3)
            Measurement zua_meas(z_zua, R_zua, id_zua, name_zua, ts_zua, mahaThresh_zua);
            ekf->insertMeasurement(zua_meas, false);

            // Zero Update Pose (ZUP): force position and attitude to stay the same (avoid drift)
            VectorX z_zup = ekf->getLastState()->x().segment(State3d::POSITION_X, 3);  // Get last estimated position
            //MatrixX R_zup = MatrixX::Identity(3, 3) * std::pow(0.3, 2.0);  // Around 1 m of 3*std
            MatrixX R_zup = ekf->getLastState()->P().block<3, 3>(State3d::POS_ERR_X, State3d::POS_ERR_X);  // Get last estimated covariance
            int id_zup = Z_POS;
            std::string name_zup = "veh_zup";
            double mahaThresh_zup = 7.81;  // chi2(0.05, 3)
            double ts_zup = currentTimestamp;
            Measurement zup_meas(z_zup, R_zup, id_zup, name_zup, ts_zup, mahaThresh_zup);
            ekf->insertMeasurement(zup_meas, false);

            // Zero Update Attitude (ZU att): force attitude to stay the same (avoid drift)
            // TODO
        } else {
            // Vehicle kinematics measurement
            double w_fl = dataset[n][13] * rpm2rads;  // [rpm] --> [rad/s]
            double w_fr = dataset[n][14] * rpm2rads;  // [rpm] --> [rad/s]
            double w_rl = dataset[n][15] * rpm2rads;  // [rpm] --> [rad/s]
            double w_rr = dataset[n][16] * rpm2rads;  // [rpm] --> [rad/s]

            double delta_f = dataset[n][17];  // Front steering angle [rad]
            double delta_r = dataset[n][18];  // Rear steering angle [rad]

            double w_mean_rads = 0.25 * (w_fl + w_fr + w_rl + w_rr);  // Mean wheel speed [rad/s]
            double vx_mean_ms = wheelRadiusM * w_mean_rads;  // [m/s]
            double num = dFrontAxleToCog * std::tan(delta_r) + dRearAxleToCog * std::tan(delta_f);
            double beta = std::atan(num / wheelBaseM);  // Sideslip angle [rad]
            double zvx = vx_mean_ms * std::cos(beta);
            double zvy = vx_mean_ms * std::sin(beta);
            double zwz = ((vx_mean_ms * std::cos(beta)) / wheelBaseM) * (std::tan(delta_f) - std::tan(delta_r));  // [rad/s]
            if (zvy > 1.0) {
                std::cerr << "WARNING: problem zvy > 1.0 m/s!" << std::endl;
            }

            VectorX z_veh = VectorX::Zero(6);
            z_veh << zvx, zvy, 0.0, 0.0, 0.0, zwz;
            MatrixX R_veh = MatrixX::Identity(6, 6);
            R_veh(0, 0) = R_veh(1, 1) = std::pow(linVelNoiseStdXY, 2.0);
            R_veh(2, 2) = std::pow(linVelNoiseStdZ, 2.0);
            R_veh(3, 3) = R_veh(4, 4) = std::pow(angVelNoiseStdXY, 2.0);
            R_veh(5, 5) = std::pow(angVelNoiseStdZ, 2.0);
            int id_veh = Z_LIN_ANG_VEL;
            std::string name_veh = "veh_kin";
            double ts_veh = currentTimestamp;
            double mahaThresh_veh = 12.59;  // chi2(0.05, 6)
            Measurement meas_veh(z_veh, R_veh, id_veh, name_veh, ts_veh, mahaThresh_veh);
            ekf->insertMeasurement(meas_veh, false);
        }

        // GPS measurement
        double gpsLatRad = dataset[n][19] * DEG_2_RAD;
        double gpsLonRad = dataset[n][20] * DEG_2_RAD;
        double gpsAltM = dataset[n][21];
        int gpsFixType = static_cast<int>(dataset[n][22]);
        int gpsNbSat = static_cast<int>(dataset[n][23]);
        
        if (gpsFixType >= 3 && gpsNbSat >= 10) {
            //std::cout << "GPS: fix type=" << gpsFixType << ", nSats=" << gpsNbSat << std::endl;
            llaToEcef(gpsLatRad, gpsLonRad, gpsAltM, &ecefX, &ecefY, &ecefZ);
            ecefToEnu(ecefX, ecefY, ecefZ, gpsLat0Rad, gpsLon0Rad, gpsAlt0M, &east, &north, &up);

            VectorX z_gps = VectorX::Zero(3);
            MatrixX R_gps = MatrixX::Identity(3, 3) * std::pow(3.0, 2.0);
            int id_gps = Z_POS;
            std::string name_gps = "gps_pos";
            double ts_gps = currentTimestamp;
            double mahaThresh_gps = 7.81;  // chi2(0.05, 3)
            Measurement gps_meas(z_gps, R_gps, id_gps, name_gps, ts_gps, mahaThresh_gps);
            ekf->insertMeasurement(gps_meas, false);
        }

        //std::cout << "Ekf last state: " << *ekf->getLastState() << std::endl;

        // Apply correction
        size_t nStates = ekf->states().size();
        if (nStates >= 2) {
            tStart = chr::high_resolution_clock::now();
            ekf->iterateFromIdx(nStates - 2);
            tStop = chr::high_resolution_clock::now();
            tSpanSec = chr::duration_cast<chr::duration<double>>(tStop - tStart).count();
            //std::cout << "Took " << tSpanSec << " s to apply correction" << std::endl;
        }

        // Sanity check
        bool isSane = ekf->checkSanity();
        if (!isSane) {
            std::bitset<8> bitState(ekf->stateSanityMask());
            std::bitset<8> bitCov(ekf->covSanityMask());
            std::cerr << "WARNING: ekf is not sane!, state mask=" << bitState << ", cov mask=" << bitCov << std::endl;
        }

        // Write results to file
        outFile << stateToString(ekf->getLastState()->x(), ekf->getLastState()->P());
        /*const VectorX x = ekf->getLastState()->x();
        outFile << x(0);
        for (int i = 1; i < x.size(); i++) {
            outFile << ";" << x(i);
        }
        outFile << "\n";*/
        
        previousTimestamp = currentTimestamp;

        n++;  // Go to next dataset row
    }  // end while

    // Save estimation filter results in a file
    outFile.close();

    delete ekf;
}  // test_ekf3

void test_ekf2(const std::vector<std::vector<double>>& dataset, const std::string& resultFilename) {
    // Output file
    std::ofstream outFile(resultFilename);

    // Create filter parameters
    FilterParameters filterParams;
    createFilterParameters(&filterParams);

    // Measurement noise covariance/std
    // --- RT2000
    const double acceleroNoiseStd = 4.0e-2; // [m/s2]
    const double gyroNoiseStd = 3.0e-3; // [rad/s]
    // --- PX4
    //const double acceleroNoiseStd = 2.0e-2; // [m/s2]
    //const double gyroNoiseStd = 2.0e-4; // [rad/s]

    const double linVelNoiseStdXY = 0.3; // [m/s]
    const double angVelNoiseStdXY = 0.05; // [rad/s]
    const double angVelNoiseStdZ = 0.03; // [rad/s]

    // Create estimation filter
    Ekf2d* ekf = new Ekf2d(filterParams);

    const int N_STATES = State2d::N_STATES;
    const int N_ERROR_STATES = State2d::N_ERROR_STATES;

    // Initialize first state of the estimation filter
    VectorX x0 = VectorX::Zero(N_STATES);
    x0(State2d::POSITION_X) = 0.0;  // 2D position
    x0(State2d::POSITION_Y) = 0.0;  // 2D position
    x0(State2d::ATTITUDE_Z) = 0.0; // Heading/yaw (psi)
    x0(State2d::LIN_VEL_X) = 0.0;  // 2D linear velocity
    x0(State2d::LIN_VEL_Y) = 0.0;  // 2D linear velocity
    x0(State2d::ANG_VEL_Z) = 0.0;  // Yaw rate
    x0(State2d::LIN_ACC_X) = 0.0;  // 2D linear acceleration
    x0(State2d::LIN_ACC_Y) = 0.0;  // 2D linear acceleration
    x0(State2d::ACC_BIAS_X) = 0.0;  // Accelero bias
    x0(State2d::ACC_BIAS_Y) = 0.0;  // Accelero bias
    x0(State2d::GYR_BIAS_Z) = 0.0;  // Gyro bias

    // Initialize covariance matrix
    MatrixX P0 = MatrixX::Identity(N_ERROR_STATES, N_ERROR_STATES);
    P0.block<2, 2>(State2d::POS_ERR_X, State2d::POS_ERR_X) = std::pow(0.3, 2.0) * Matrix2::Identity();
    P0(State2d::ATT_ERR_Z, State2d::ATT_ERR_Z) = std::pow(0.03, 2.0);
    P0.block<2, 2>(State2d::LIN_VEL_ERR_X, State2d::LIN_VEL_ERR_X) = std::pow(0.1, 2.0) * Matrix2::Identity();
    P0(State2d::ANG_VEL_ERR_Z, State2d::ANG_VEL_ERR_Z) = std::pow(0.01, 2.0);
    P0.block<2, 2>(State2d::LIN_ACC_ERR_X, State2d::LIN_ACC_ERR_X) = std::pow(0.01, 2.0) * Matrix2::Identity();
    P0.block<2, 2>(State2d::ACC_BIAS_ERR_X, State2d::ACC_BIAS_ERR_X) = std::pow(0.15, 2.0) * Matrix2::Identity();
    P0(State2d::GYR_BIAS_ERR_Z, State2d::GYR_BIAS_ERR_Z) = std::pow(0.05, 2.0);

    const double t0 = dataset[0][0];  // [s]

    // Initialize filter
    ekf->init(x0, P0, t0);
    outFile << stateToString(ekf->getLastState()->x(), ekf->getLastState()->P());  // Add to result

    // GPS initialization
    double gpsLat0Rad = dataset[0][19] * DEG_2_RAD;
    double gpsLon0Rad = dataset[0][20] * DEG_2_RAD;
    double gpsAlt0M = dataset[0][21];
    double ecefX = 0.0, ecefY = 0.0, ecefZ = 0.0;
    double east = 0.0, north = 0.0, up = 0.0;
    llaToEcef(gpsLat0Rad, gpsLon0Rad, gpsAlt0M, &ecefX, &ecefY, &ecefZ);
    ecefToEnu(ecefX, ecefY, ecefZ, gpsLat0Rad, gpsLon0Rad, gpsAlt0M, &east, &north, &up);

    // Initialize useful variables
    double previousTimestamp = t0;  // [s]
    double currentTimestamp = t0;  // [s]
    const double timestampThresh = 0.01;  // [s]
    bool vehicleAtRest = false;  // If the vehicle is stopped and not moving

    const double wheelRadiusM = 0.3575;  // [m]
    const double dFrontAxleToCog = 1.135;  // [m]
    const double dRearAxleToCog = 1.115;  // [m]
    const double wheelBaseM = dFrontAxleToCog + dRearAxleToCog;  // [m]
    const double rpm2rads = M_PI / 30.0;

    // Timing variables
    chr::high_resolution_clock::time_point tStart, tStop;
    double tSpanSec = 0.0;  // [s]

    std::cout << "Main loop..." << std::endl;

    // Main temporal loop over the dataset
    //bool isRunning = true;
    const size_t nRows = dataset.size();
    size_t n = 1;
    //while (isRunning) {
    while (n < nRows) {
        currentTimestamp = dataset[n][0];  // [s]
        if (0 == n % 100)
            std::cout << "Iteration " << n << ": time=" << currentTimestamp << " s" << std::endl;

        // Create new state and add it to the filter
        double deltaTimestamp = currentTimestamp - ekf->getLastState()->timestamp();
        //std::cout << "Last state timestamp=" << ekf->getLastState()->timestamp() << " s" << std::endl;
        if (deltaTimestamp >= timestampThresh) {
            //std::cerr << "Creating a new state..." << std::endl;
            State2d state(currentTimestamp);
            //tStart = chr::high_resolution_clock::now();
            ekf->insertState(state);
            //tStop = chr::high_resolution_clock::now();
            //tSpanSec = chr::duration_cast<chr::duration<double>>(tStop - tStart).count();
            //std::cout << "Took " << tSpanSec << " s to insert a new state" << std::endl;
        }

        // Check for new measurements

        // IMU measurement
        double zax = dataset[n][1], zay = dataset[n][2];  // RT2000 accelero
        double zgz = dataset[n][6];  // RT2000 gyro

        double zax2 = dataset[n][7], zay2 = dataset[n][8];  // PX4 accelero
        double zgz2 = dataset[n][12];  // PX4 gyro

        VectorX z_imu = VectorX::Zero(3);
        //z_imu << zax, zay, zgz;
        z_imu << zax2, zay2, zgz2;
        MatrixX R_imu = MatrixX::Identity(3, 3);
        R_imu(0, 0) = R_imu(1, 1) = std::pow(acceleroNoiseStd, 2.0);
        R_imu(2, 2) = std::pow(gyroNoiseStd, 2.0);
        int id_imu = Z_IMU;
        std::string name_imu = "imu";
        double ts_imu = currentTimestamp;
        double mahaThresh_imu = 7.81;  // chi2(0.05, 3)
        Measurement meas_imu(z_imu, R_imu, id_imu, name_imu, ts_imu, mahaThresh_imu);
        ekf->insertMeasurement(meas_imu, false);

        // Vehicle measurements
        int transmissionMode = static_cast<int>(dataset[n][24]);  // 0: sna, 1: park, 2: reverse, 3: neutral, 4: forward
        double vehicleSpeedKmh = dataset[n][25];  // [km/h]

        vehicleAtRest = (((1 == transmissionMode) || (3 == transmissionMode)) && (std::abs(vehicleSpeedKmh) <= 0.1)) ? true : false;

        if (vehicleAtRest) {  // Vehicle at rest: ZUP, ZUV, ZUA
            // Zero Update Velocities (ZUV): force linear and angular velocities to zero
            VectorX z_zuv = VectorX::Zero(3);
            MatrixX R_zuv = MatrixX::Identity(3, 3);
            R_zuv(0, 0) = R_zuv(1, 1) = std::pow(0.01, 2.0);  // Around 0.1 km/h 3*std
            R_zuv(2, 2) = std::pow(0.001, 2.0);  // Around 0.2 deg/s 3*std
            int id_zuv = Z_LIN_ANG_VEL;
            std::string name_zuv = "veh_zuv";
            double ts_zuv = currentTimestamp;
            double mahaThresh_zuv = 7.81;  // chi2(0.05, 3)
            Measurement meas_zuv(z_zuv, R_zuv, id_zuv, name_zuv, ts_zuv, mahaThresh_zuv);
            ekf->insertMeasurement(meas_zuv, false);

            // Zero Update Acceleration (ZUA): force linear acceleration to zero
            VectorX z_zua = VectorX::Zero(2);
            MatrixX R_zua = MatrixX::Identity(2, 2) * std::pow(0.0033, 2.0);  // Around 0.01 m/s2 3*std
            int id_zua = Z_LIN_ACC;
            std::string name_zua = "veh_zua";
            double ts_zua = currentTimestamp;
            double mahaThresh_zua = 5.99;  // chi2(0.05, 2)
            Measurement zua_meas(z_zua, R_zua, id_zua, name_zua, ts_zua, mahaThresh_zua);
            ekf->insertMeasurement(zua_meas, false);

            // Zero Update Pose (ZUP): force position and attitude to stay the same (avoid drift)
            VectorX z_zup = ekf->getLastState()->x().segment(State2d::POSITION_X, 2);  // Get last estimated position
            //MatrixX R_zup = MatrixX::Identity(2, 2) * std::pow(0.3, 2.0);  // Around 1 m of 3*std
            MatrixX R_zup = ekf->getLastState()->P().block<2, 2>(State2d::POS_ERR_X, State2d::POS_ERR_X);  // Get last estimated covariance
            int id_zup = Z_POS;
            std::string name_zup = "veh_zup";
            double mahaThresh_zup = 5.99;  // chi2(0.05, 2)
            double ts_zup = currentTimestamp;
            Measurement zup_meas(z_zup, R_zup, id_zup, name_zup, ts_zup, mahaThresh_zup);
            ekf->insertMeasurement(zup_meas, false);

            // Zero Update Attitude (ZU att): force attitude to stay the same (avoid drift)
            // TODO
        } else {
            // Vehicle kinematics measurement
            double w_fl = dataset[n][13] * rpm2rads;  // [rpm] --> [rad/s]
            double w_fr = dataset[n][14] * rpm2rads;  // [rpm] --> [rad/s]
            double w_rl = dataset[n][15] * rpm2rads;  // [rpm] --> [rad/s]
            double w_rr = dataset[n][16] * rpm2rads;  // [rpm] --> [rad/s]

            double delta_f = dataset[n][17];  // Front steering angle [rad]
            double delta_r = dataset[n][18];  // Rear steering angle [rad]

            double w_mean_rads = 0.25 * (w_fl + w_fr + w_rl + w_rr);  // Mean wheel speed [rad/s]
            double vx_mean_ms = wheelRadiusM * w_mean_rads;  // [m/s]
            double num = dFrontAxleToCog * std::tan(delta_r) + dRearAxleToCog * std::tan(delta_f);
            double beta = std::atan(num / wheelBaseM);  // Sideslip angle [rad]
            double zvx = vx_mean_ms * std::cos(beta);
            double zvy = vx_mean_ms * std::sin(beta);
            double zwz = ((vx_mean_ms * std::cos(beta)) / wheelBaseM) * (std::tan(delta_f) - std::tan(delta_r));  // [rad/s]
            if (zvy > 1.0) {
                std::cerr << "WARNING: problem zvy > 1.0 m/s!" << std::endl;
            }

            VectorX z_veh = VectorX::Zero(3);
            z_veh << zvx, zvy, zwz;
            MatrixX R_veh = MatrixX::Identity(3, 3);
            R_veh(0, 0) = R_veh(1, 1) = std::pow(linVelNoiseStdXY, 2.0);
            R_veh(2, 2) = std::pow(angVelNoiseStdZ, 2.0);
            int id_veh = Z_LIN_ANG_VEL;
            std::string name_veh = "veh_kin";
            double ts_veh = currentTimestamp;
            double mahaThresh_veh = 7.81;  // chi2(0.05, 3)
            Measurement meas_veh(z_veh, R_veh, id_veh, name_veh, ts_veh, mahaThresh_veh);
            ekf->insertMeasurement(meas_veh, false);
        }

        // GPS measurement
        double gpsLatRad = dataset[n][19] * DEG_2_RAD;
        double gpsLonRad = dataset[n][20] * DEG_2_RAD;
        double gpsAltM = dataset[n][21];
        int gpsFixType = static_cast<int>(dataset[n][22]);
        int gpsNbSat = static_cast<int>(dataset[n][23]);
        
        if (gpsFixType >= 3 && gpsNbSat >= 10) {
            //std::cout << "GPS: fix type=" << gpsFixType << ", nSats=" << gpsNbSat << std::endl;
            llaToEcef(gpsLatRad, gpsLonRad, gpsAltM, &ecefX, &ecefY, &ecefZ);
            ecefToEnu(ecefX, ecefY, ecefZ, gpsLat0Rad, gpsLon0Rad, gpsAlt0M, &east, &north, &up);

            VectorX z_gps = VectorX::Zero(2);
            MatrixX R_gps = MatrixX::Identity(2, 2) * std::pow(3.0, 2.0);
            int id_gps = Z_POS;
            std::string name_gps = "gps_pos";
            double ts_gps = currentTimestamp;
            double mahaThresh_gps = 5.99;  // chi2(0.05, 2)
            Measurement gps_meas(z_gps, R_gps, id_gps, name_gps, ts_gps, mahaThresh_gps);
            ekf->insertMeasurement(gps_meas, false);
        }

        //std::cout << "Ekf last state: " << *ekf->getLastState() << std::endl;

        // Apply correction
        size_t nStates = ekf->states().size();
        if (nStates >= 2) {
            tStart = chr::high_resolution_clock::now();
            ekf->iterateFromIdx(nStates - 2);
            tStop = chr::high_resolution_clock::now();
            tSpanSec = chr::duration_cast<chr::duration<double>>(tStop - tStart).count();
            //std::cout << "Took " << tSpanSec << " s to apply correction" << std::endl;
        }

        // Sanity check
        bool isSane = ekf->checkSanity();
        if (!isSane) {
            std::bitset<8> bitState(ekf->stateSanityMask());
            std::bitset<8> bitCov(ekf->covSanityMask());
            std::cerr << "WARNING: ekf is not sane!, state mask=" << bitState << ", cov mask=" << bitCov << std::endl;
        }

        // Write results to file
        outFile << stateToString(ekf->getLastState()->x(), ekf->getLastState()->P());
        /*const VectorX x = ekf->getLastState()->x();
        outFile << x(0);
        for (int i = 1; i < x.size(); i++) {
            outFile << ";" << x(i);
        }
        outFile << "\n";*/
        
        previousTimestamp = currentTimestamp;

        n++;  // Go to next dataset row
    }  // end while

    // Save estimation filter results in a file
    outFile.close();

    delete ekf;
}  // test_ekf2

int main(int argc, char** argv) {
    // Parse arguments if any
    std::string datasetFilename = "../Data/ekfDataset.csv";
    std::string resultFilename = "../Data/ekfResults.csv";
    if (argc > 1) {
        datasetFilename = std::string(argv[1]);
        std::cout << "Given dataset filename is: " << datasetFilename << std::endl;
    }
    if (argc > 2) {
        resultFilename = std::string(argv[2]);
        std::cout << "Given result filename is: " << resultFilename << std::endl;
    }

    // Load dataset
    std::vector<std::vector<double>> dataset;
    if (!loadDataset(datasetFilename, dataset)) {
        std::cerr << "Failed to load dataset from file: " << datasetFilename << std::endl;
        return -1;
    }

    std::string firstRow = row2str(dataset[0]);
    std::string lastRow = row2str(dataset.back());
    std::cout << "First dataset row: " << firstRow << std::endl;
    std::cout << "Last dataset row: " << lastRow << std::endl;

    //test_ekf3(dataset, resultFilename);
    test_ekf2(dataset, resultFilename);

    return 0;
}  // main