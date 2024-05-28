///
/// \file EkfLocNode.cpp
/// \brief Ekf localization ROS node with GPS, odometry and IMU measurements
/// \author Lionel GENEVE
/// \date 20/11/2018
/// \version 1.0
///

// STD
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <queue>
#include <utility>
#include <cmath>

// EIGEN
#include <Eigen/Dense>

// ROS
#include "ros/ros.h"
#include "geometry_msgs/PointStamped.h"
#include "geometry_msgs/PoseStamped.h"
#include "geometry_msgs/PoseWithCovarianceStamped.h"
#include "geometry_msgs/TwistStamped.h"
#include "geometry_msgs/TwistWithCovarianceStamped.h"
#include "sensor_msgs/Imu.h"
#include "nav_msgs/Odometry.h"
//#include "nav_msgs/Path.h"


#include "Ekf.hpp"

// Global variables
static double gpsPositionMahaThresh = 10.597; // df=3, p=0.005
static double gpsPoseMahaThresh = 12.838;
static double gpsOrientationMahaThresh = 8.;
static double odomLinAngVelMahaThresh = 10.597; //  df=2, p=0.005
static double odomLinVelMahaThresh = 100.;
static double odomAngVelMahaThresh = 100.;

std::queue<ekf::Measurement> measurementQueue;
bool isInit = false;

///
/// \fn gpsCallback
/// \brief ROS callback function for the GPS message
///
void gpsCallback(const geometry_msgs::PoseWithCovarianceStamped::ConstPtr& _msg) {
    //std::cerr << "EkfLoc: gpsCallback\n";
    double tStamp = _msg->header.stamp.toSec();
    std::string sensorName = "gps";
    Eigen::VectorXd z;
    Eigen::MatrixXd R;
    int id = -1;
    double maha = gpsPositionMahaThresh;

    double covHeading = _msg->pose.covariance[35];
    if (covHeading < 0) {
        id = ekf::MeasurementIdx_t::Z_PX_PY;
        z = Eigen::VectorXd::Zero(2);
        R = Eigen::MatrixXd::Zero(2,2);
        z << _msg->pose.pose.position.x, _msg->pose.pose.position.y;
        R << _msg->pose.covariance[0], _msg->pose.covariance[1],
            _msg->pose.covariance[6], _msg->pose.covariance[7];
        maha = gpsPositionMahaThresh;
    }
    else {
        if (!isInit) {
            id = ekf::MeasurementIdx_t::Z_PX_PY_PT;
            z = Eigen::VectorXd::Zero(3);
            R = Eigen::MatrixXd::Zero(3,3);
            double theta = 2. * std::atan2(_msg->pose.pose.orientation.z,
                    _msg->pose.pose.orientation.w);
            z << _msg->pose.pose.position.x, _msg->pose.pose.position.y, theta;
            R << _msg->pose.covariance[0], _msg->pose.covariance[1], 0.,
                _msg->pose.covariance[6], _msg->pose.covariance[7], 0.,
                0., 0., covHeading;
            maha = gpsPoseMahaThresh;
        }
        else {
            // Separate GPS orientation and position measurement
            id = ekf::MeasurementIdx_t::Z_PT;
            z = Eigen::VectorXd::Zero(1);
            R = Eigen::MatrixXd::Zero(1,1);
            z << 2. * std::atan2(_msg->pose.pose.orientation.z,
                    _msg->pose.pose.orientation.w);
            R << covHeading;
            ekf::Measurement measT(z, R, id, sensorName, tStamp, gpsOrientationMahaThresh);
            measurementQueue.push(measT);

            id = ekf::MeasurementIdx_t::Z_PX_PY;
            z = Eigen::VectorXd::Zero(2);
            R = Eigen::MatrixXd::Zero(2,2);
            z << _msg->pose.pose.position.x, _msg->pose.pose.position.y;
            R << _msg->pose.covariance[0], _msg->pose.covariance[1],
                _msg->pose.covariance[6], _msg->pose.covariance[7];
            maha = gpsPositionMahaThresh;
        }
    }
    
    ekf::Measurement meas(z, R, id, sensorName, tStamp, maha);

    measurementQueue.push(meas);
}

///
/// \fn odomCallback
/// \brief ROS callback function for the odometry message
///
void odomCallback(const nav_msgs::Odometry::ConstPtr& _msg) {
    //std::cerr << "EkfLoc: odomCallback\n";
    double tStamp = _msg->header.stamp.toSec();
    std::string sensorName = "odom";
    
    //Eigen::VectorXd z = Eigen::VectorXd::Zero(2);
    //Eigen::MatrixXd R = Eigen::MatrixXd::Zero(2,2);
    //int id = ekf::MeasurementIdx_t::Z_VX_WZ;
    //z << _msg->twist.twist.linear.x, _msg->twist.twist.angular.z;
    //R << _msg->twist.covariance[0], 0.,
        //0., _msg->twist.covariance[35];
    //ekf::Measurement meas(z, R, id, sensorName, tStamp, odomLinAngVelMahaThresh);
    //measurementQueue.push(meas);

    // Separate linear and angular odometry measurement
    Eigen::VectorXd z = Eigen::VectorXd::Zero(1);
    Eigen::MatrixXd R = Eigen::MatrixXd::Zero(1,1);
    int id = ekf::MeasurementIdx_t::Z_VX;
    z << _msg->twist.twist.linear.x;
    R << _msg->twist.covariance[0];
    ekf::Measurement measV(z, R, id, sensorName, tStamp, odomLinVelMahaThresh);
    measurementQueue.push(measV);

    id = ekf::MeasurementIdx_t::Z_WZ;
    z << _msg->twist.twist.angular.z;
    R << _msg->twist.covariance[35];
    ekf::Measurement measW(z, R, id, sensorName, tStamp, odomAngVelMahaThresh);
    measurementQueue.push(measW);
}

///
/// \fn imuCallback
/// \brief ROS callback function for the IMU message
///
void imuCallback(const sensor_msgs::Imu::ConstPtr& /*_msg*/) {
    //std::cerr << "EkfLoc: imuCallback\n";
}

///
/// \fn main
/// \brief main function
///
int main(int argc, char** argv) {
    ros::init(argc, argv, "ekf_loc_node");
    // ROS handlers
    ros::NodeHandle localNh("~"); // Local ROS node handler
    ros::NodeHandle globalNh(""); // Global ROS node handler

    // ROS parameters
    std::string odomFrame = "map";
    std::string robotFrame = "robot";
    int nIters = 1;
    double temporalWindow = 10.;  // Temporal window of the buffer of states [s]

    localNh.param<std::string>("odomFrame", odomFrame, "odom");
    localNh.param<std::string>("robotFrame", robotFrame, "robot");
    localNh.param<int>("nIters", nIters, 1);
    localNh.param<double>("temporalWindow", temporalWindow, 10.);
    localNh.param<double>("gpsPositionMahaThresh", gpsPositionMahaThresh, 100.);
    localNh.param<double>("gpsPoseMahaThresh", gpsPoseMahaThresh, 100.);
    localNh.param<double>("gpsOrientationMahaThresh", gpsOrientationMahaThresh, 100.);
    localNh.param<double>("odomLinAngVelMahaThresh", odomLinAngVelMahaThresh, 100.);
    localNh.param<double>("odomLinVelMahaThresh", odomLinVelMahaThresh, 100.);
    localNh.param<double>("odomAngVelMahaThresh", odomAngVelMahaThresh, 100.);

    std::cerr << "EkfLoc parameters:\n";
    std::cerr << "\t-odomFrame=" << odomFrame << std::endl;
    std::cerr << "\t-robotFrame=" << robotFrame << std::endl;
    std::cerr << "\t-nIters=" << nIters << std::endl;
    std::cerr << "\t-temporalWindow=" << temporalWindow << std::endl;
    std::cerr << "\t-gpsPositionMahaThresh=" << gpsPositionMahaThresh << std::endl;
    std::cerr << "\t-gpsPoseMahaThresh=" << gpsPoseMahaThresh << std::endl;
    std::cerr << "\t-gpsOrientationMahaThresh=" << gpsOrientationMahaThresh << std::endl;
    std::cerr << "\t-odomLinAngVelMahaThresh=" << odomLinAngVelMahaThresh << std::endl;
    std::cerr << "\t-odomLinVelMahaThresh=" << odomLinVelMahaThresh << std::endl;
    std::cerr << "\t-odomAngVelMahaThresh=" << odomAngVelMahaThresh << std::endl;

    // ROS subscribers
    ros::Subscriber gpsSub =
        globalNh.subscribe<geometry_msgs::PoseWithCovarianceStamped>("gps", 1,
        &gpsCallback);
    ros::Subscriber odomSub =
        globalNh.subscribe<nav_msgs::Odometry>("odom", 1, &odomCallback);
    ros::Subscriber imuSub =
        globalNh.subscribe<sensor_msgs::Imu>("imu", 10, &imuCallback);

    // ROS publishers
    ros::Publisher odomPub =
        globalNh.advertise<nav_msgs::Odometry>("ekf/odom", 1);
    ros::Publisher posePub =
        globalNh.advertise<geometry_msgs::PoseWithCovarianceStamped>("ekf/pose", 1);

    // ROS messages
    nav_msgs::Odometry odomMsg;
    odomMsg.header.frame_id = odomFrame;
    odomMsg.header.seq = 0;
    odomMsg.child_frame_id = robotFrame;
    for (int i=0;i<36;i++) {
        odomMsg.pose.covariance[i] = 0.;
        odomMsg.twist.covariance[i] = 0.;
    }

    ekf::Ekf ekf(nIters, temporalWindow);

    std::vector<double> vs, ws;
    std::vector<ekf::Measurement> odomMeas;
    bool publish = false;
    isInit = false;
    ros::Rate loopRate(500); // 500 hz

    while (ros::ok()) {
        if (!measurementQueue.empty()) { // If new measurement
            ekf::Measurement meas = measurementQueue.front(); // Retrieve new measurement
            measurementQueue.pop(); // Delete from the queue
            //std::cerr << "EkfLoc: pop measurement: " << meas << std::endl;

            if (ekf.isInit()) {
                //std::cerr << "EkfLoc: inserting measurement in the filter\n";
                ekf.insertMeasurement(meas);
                publish = true;
            }
            else {
                //std::cerr << "Ekf loc filter not initialized\n";
                // Try to initialize the filter
                if (meas.id() == ekf::MeasurementIdx_t::Z_PX_PY_PT) {
                    std::cerr << "EkfLoc: received 1st pose measurement from GPS\n";
                    double tStamp = meas.timestamp();
                    int nStates = ekf::StateIdx_t::N_STATES;
                    Eigen::VectorXd x0 = Eigen::VectorXd::Zero(nStates);
                    Eigen::MatrixXd P0 = Eigen::MatrixXd::Zero(nStates,nStates);
                    double v = 0., w = 0.;
                    if (vs.size() > 0.)
                        v = vs.back();
                    if (ws.size() > 0.)
                        w = ws.back();
                    x0.block<3,1>(0,0) = meas.z();
                    x0.block<2,1>(3,0) << v, w;
                    P0.block<3,3>(0,0) = meas.R();
                    P0(2,2) *= 2.;
                    P0(3,3) = std::pow(.5 / 3., 2.);
                    P0(4,4) = std::pow((1. * M_PI / 180.) / 3., 2.);
                    std::cerr << "EkfLoc: x0=" << x0.transpose() << std::endl;
                    std::cerr << "EkfLoc: P0=\n" << P0 << std::endl;
                    ekf.init(x0, P0, tStamp);
                    if (odomMeas.size() > 0) {
                        size_t idx = odomMeas.size() - 1;
                        while ((idx > 0) && odomMeas[idx].timestamp() > tStamp) {
                            ekf.insertMeasurement(odomMeas[idx]);
                            odomMeas.erase(odomMeas.begin()+idx);
                            idx--;
                        }
                        odomMeas.clear();
                        vs.clear();
                        ws.clear();
                    }
                    publish = true;
                    isInit = true;
                }
                else if (meas.id() == ekf::MeasurementIdx_t::Z_VX_WZ) {
                    //std::cerr << "Accumulating odom measurement\n";
                    vs.push_back(meas.z()(0));
                    ws.push_back(meas.z()(1));
                    odomMeas.push_back(meas);
                }
                else if (meas.id() == ekf::MeasurementIdx_t::Z_VX) {
                    vs.push_back(meas.z()(0));
                    odomMeas.push_back(meas);
                }
                else if (meas.id() == ekf::MeasurementIdx_t::Z_WZ) {
                    ws.push_back(meas.z()(0));
                    odomMeas.push_back(meas);
                }
            }
        }

        if (publish) {
            //std::cerr << "EkfLoc: Publishing ekf state\n";
            publish = false;
            ekf::State* state = ekf.getLastState();
            if (state != nullptr) {
                odomMsg.header.stamp = ros::Time(state->timestamp());
                odomMsg.header.seq++;
                // Position
                odomMsg.pose.pose.position.x = state->x()(ekf::StateIdx_t::PX);
                odomMsg.pose.pose.position.y = state->x()(ekf::StateIdx_t::PY);
                //odomMsg.pose.pose.position.z = 0.;
                // Orientation
                double theta = state->x()(ekf::StateIdx_t::PT);
                odomMsg.pose.pose.orientation.w = std::cos(.5 * theta);
                //odomMsg.pose.pose.orientation.x = 0.;
                //odomMsg.pose.pose.orientation.y = 0.;
                odomMsg.pose.pose.orientation.z = std::sin(.5 * theta);
                // Linear velocity
                odomMsg.twist.twist.linear.x = state->x()(ekf::StateIdx_t::VX);
                //odomMsg.twist.twist.linear.y = 0.;
                //odomMsg.twist.twist.linear.z = 0.;
                // Angular velocity
                //odomMsg.twist.twist.angular.x = 0.;
                //odomMsg.twist.twist.angular.y = 0.;
                odomMsg.twist.twist.angular.z = state->x()(ekf::StateIdx_t::WZ);
                // Pose covariance
                odomMsg.pose.covariance[0] = state->P()(ekf::PX,ekf::PX);
                odomMsg.pose.covariance[1] = state->P()(ekf::PX,ekf::PY);
                odomMsg.pose.covariance[5] = state->P()(ekf::PX,ekf::PT);
                odomMsg.pose.covariance[6] = state->P()(ekf::PY,ekf::PX);
                odomMsg.pose.covariance[7] = state->P()(ekf::PY,ekf::PY);
                odomMsg.pose.covariance[30] = state->P()(ekf::PT,ekf::PX);
                odomMsg.pose.covariance[35] = state->P()(ekf::PT,ekf::PT);
                // Velocity covariance
                odomMsg.twist.covariance[0] = state->P()(ekf::VX,ekf::VX);
                odomMsg.twist.covariance[35] = state->P()(ekf::WZ,ekf::WZ);
                
                odomPub.publish(odomMsg);

                //std::cerr << "EkfLoc: P=\n" << state->P() << std::endl;
            }
        }

        loopRate.sleep();
        ros::spinOnce();
    }

    ros::waitForShutdown();

    return 0;
}
