function particles = fastslam_propagate_2d_SLAM(particles, v_m, omega_m, dt, sigma_v, sigma_omega);
% fastslam_propagate_2d_SLAM: sampling a new pose, as in FastSLAM v1.0

particleNum = length(particles);

% generate noise
v = zeros(particleNum,1);
omega = zeros(particleNum,1);

v = v_m - sigma_v*randn(particleNum,1);
omega = omega_m - sigma_omega*randn(particleNum,1);

for i=1:particleNum
    pose = particles(i).Xr;
    particles(i).Xr(1) = pose(1) + v(i)*dt*cos(pose(3));
    particles(i).Xr(2) = pose(2) + v(i)*dt*sin(pose(3));
    particles(i).Xr(3) = pose(3) + omega(i)*dt;
end