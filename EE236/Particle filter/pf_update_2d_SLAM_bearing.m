function particles = pf_update_2d_SLAM_bearing(particles, z, sigma_theta, X_L)
% Update the weights of particles based on measurements

particleNum = length(particles);

for i=1:particleNum
    pose = particles(i).Xr;
    C = [cos(pose(3)) -sin(pose(3));
        sin(pose(3)) cos(pose(3))];
    
    R_P_L = C'*(X_L - pose(1:2));
    
    % the estimated measurement
    % the range and bearing:
	[z_hat, range] = cart2pol(R_P_L(1),R_P_L(2));
    % measurement innovation
    r = z - z_hat; 

    % Calculate the importance weight
    % Equation (62), (63). Denominator is not needed?
    gauss = exp(-0.5*r^2/sigma_theta^2)/sqrt(2*pi*sigma_theta^2);
    particles(i).weight = particles(i).weight*gauss;
end

% Weights will be normalized in the resampling function
        
end