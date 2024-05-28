# StateEstimationFilter

## Description
State estimation filter using an error-state extended Kalman filter (ES-EKF).
The filter maintains a temporal sliding window of past states to integrate delayed measurements.

## State vector
The state vector **x** being estimated is composed by:
- the 3D position **p** [m] in local ENU frame,
- the 3D attitude (quaternion) **q** from body frame to local ENU frame,
- the 3D linear velocity **v** [m/s] in body frame,
- the 3D angular velocity **w** [rad/s] in body frame,
- the 3D linear acceleration **a** [m/s2] in body frame,
- the accelerometer bias **b_a** [m/s2] in sensor frame,
- the gyroscope bias **b_g** [rad/s] in sensor frame,
- optional: the magnetometer bias **b_m** [uT] in sensor frame (not implemented)

## Prediction
The filter prediction/propagation step is based on a constant linear acceleration and angular velocity model.

## Measurements
The measurements used to correct the estimated states are:
- IMU: accelerometer and gyroscope
- Magnetometer (optional)
- Vehicle kinematics constraints on the linear and angular velocities
- Zero update velocity/position/acceleration when the vehicle is at rest
- GPS (position, velocity, yaw)
- others (lidar, camera...)

## References
List of references:
- Quaternion kinematics for the error-state Kalman, J. Sola
- Optimal state estimation, D. Simon
- Probabilistic Robotics, S. Thrun, W. Burgard, D. Fox

## License
This project is licensed under the Apache 2.0 license.
