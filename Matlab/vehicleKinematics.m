function [v ,w] = vehicleKinematics(wheel_speed_rpm, steer_angle_rad, ...
    wheel_radius, dist_cog_front_axle, dist_cog_rear_axle)
%vehicleKinematics: compute the vehicle velocities using the vehicle kinematics model
% Assumptions: four-wheel drive vehicle with two steering (front and rear axle)
% Args:
% - wheel_speed_rpm (vector [4x1]): wheel speeds in [rpm] (order: FL, FR, RL, RR)
% - steer_angle_rad (vector [2x1]): steering/wheel angles in [rad] (order: front, rear)
% - wheel_radius (double): tyre/wheel radius in [m] (assumed same for all wheels)
% - dist_cog_front_axle (double): distance between vehicle CoG and front axle in [m]
% - dist_cog_rear_axle (double): distance between vehicle CoG and rear axle in [m]
% Return:
% - v (vector [3x1]): linear velocity vector [m/s] v=(v_x, v_y, v_z)
% - w (vector [3x1]): angular velocity vector [rad/s] w=(w_x, w_y, w_z)

rpm2rads = pi / 30.0;

% Wheel speed in [m/s]
v_fl = wheel_speed_rpm(1) * rpm2rads * wheel_radius; % [m/s]
v_fr = wheel_speed_rpm(2) * rpm2rads * wheel_radius; % [m/s]
v_rl = wheel_speed_rpm(3) * rpm2rads * wheel_radius; % [m/s]
v_rr = wheel_speed_rpm(4) * rpm2rads * wheel_radius; % [m/s]

delta_f = steer_angle_rad(1); % [rad]
delta_r = steer_angle_rad(2); % [rad]

% Vehicle linear speed at CoG [m/s]
vx = 0.25 * (v_fl .* cos(delta_f) + v_fr .* cos(delta_f) ...
    + v_rl .* cos(delta_r) + v_rr .* cos(delta_r));
vy = 0.25 * (v_fl .* sin(delta_f) + v_fr .* sin(delta_f) ...
    + v_rl .* sin(delta_r) + v_rr .* sin(delta_r));
vz = 0.0;

% Vehicle angular speed at CoG [rad/s]
wx = 0.0;
wy = 0.0;
wheel_base = dist_cog_front_axle + dist_cog_rear_axle;
num = dist_cog_front_axle * tan(delta_r) + dist_cog_rear_axle * tan(delta_f);
beta = atan(num / wheel_base);
wz = ((vx .* cos(beta)) ./ wheel_base) .* (tan(delta_f) - tan(delta_r));

% Set output
v = [vx; vy; vz];
w = [wx; wy; wz];

end
