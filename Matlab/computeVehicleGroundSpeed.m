function v_ms = computeVehicleGroundSpeed(wheel_speed_rpm, wheel_radius)
% computeVehicleGroundSpeed: compute vehicle absolute ground speed from wheel speeds
% Args:
% - wheel_speed_rpm (vector [4x1]): speed of each wheel [rpm] (order: FL, FR, RL, RR)
% - wheel_radius (double): wheel/tyre rolling radius [m] (assumed same for all wheels)
% Return:
% - v_ms (double): vehicle ground speed [m/s]

rpm2rads = pi / 30.0;
v_ms = mean(wheel_speed_rpm) * rpm2rads * wheel_radius;

end