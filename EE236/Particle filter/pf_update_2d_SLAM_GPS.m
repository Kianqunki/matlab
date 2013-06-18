function particles = pf_update_2d_SLAM_GPS(particles, z, R_GPS)
% Update the weights of particles based on measurements

particleNum = length(particles);

for i=1:particleNum
    pose = particles(i).Xr;

    % the estimated GPS measurement
    z_hat = pose(1:2,:);
    % measurement innovation
    r = z - z_hat; 

    % Calculate the importance weight
    % Equation (62), (63). Denominator is not needed?
    gauss = exp(-0.5*r'*inv(R_GPS)*r)/(2*pi*sqrt(det(R_GPS)));
    particles(i).weight = particles(i).weight*gauss;
end

% Weights will be normalized in the resampling function
        
end