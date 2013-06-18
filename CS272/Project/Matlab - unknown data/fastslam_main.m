clear all;
close all;

%set random seed
seedRandom = 1;
rand('state',seedRandom);
randn('state',seedRandom);

% setting up the simulation

% Parameters for process model and measurement model
% robot's velocity v
vel = 0.3; % m/sec
% robot's rotational velocity omega
omega = 0.015; % rad/sec
% time length of simulation, aka number of timesteps
timesteps = 1000;
% sampling interval for sensors, aka delta_t
dt=1;
% std. deviation for odometry measurements:
sigma_vel = 0.02; % m/sec
sigma_omega = 0.002; % rad/sec

% covariance matrix of landmark measurements
%R = 0.01*eye(2); % m^2
% In FastSLAM 1.0, if measurement error is small, it diverges
% in thesis page 74
R = 0.5*[1 0.5; 0.5 1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters for simulated world
% landmarks on an L_N x L_N grid 
L_N= 9;

% range of laser scanner
%(not all landmarks will be within the robot's sensing range)
max_range = 10; %m
min_range = 0.5; %m

% generate the simulated world
[X, v_m, omega_m, z, X_L] = simulated_world_generator(vel, omega, timesteps, dt, sigma_vel, sigma_omega, R, min_range, max_range, L_N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Paramters for FastSLAM
% number of particles M
numParticle = 10;
% minimum effective particle
Nmin = numParticle*0.6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the filter:

%initialize the particles
% particles = repmat(struct('weight',1/numParticle,'Xr', [0;0;0],'Xl', [],'Pl',[]),numParticle, 1);
particles = repmat(struct('weight',1/numParticle,...
    'Xr', X(:,1),...  % robot pose: position and orientation
    'Xl', [],...       % matrix of landmark positions, 2xn
    'Pl', [],...        % matrix of landmark covariances
    'dict', []),...        % DEBUG: true index of landmarks, for debugging
    numParticle, 1);

% we also maintain a vector with the descriptors of the landmarks that we
% have added to the state vector:
descriptors = [];  % initially we don't have any landmarks

% save the robot positions for plotting trajectory:
Xr_sav = zeros(size(X,1),timesteps);
Xr_sav(:,1) = X(:,1);

% starting SLAM
for k=1:timesteps-1
    % Read in the measurements from the RWS
    
    % propagation
    particles = fastslam_propagate_2d_SLAM(particles, v_m(k), omega_m(k), dt, sigma_vel, sigma_omega);
    
    % update, with unknown data associtation
    [particles, descriptors] = fastslam_update_associate_2d_SLAM(particles, z(k+1), R, descriptors,X_L);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Importance resampling
    % section 3.3.4 in the thesis
    particles = fastslam_resample( particles, Nmin );
     	
	% save these quantities for plotting
    weight_all = [particles.weight];
    [w,idx] = max(weight_all);
    Xr_sav(:,k+1) = particles(idx).Xr(:);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotting:

% compute the estimation error for the robot:
err = X - Xr_sav;
% make sure the orientation error is in the first cycle:
err(3,:) = atan2( sin(err(3,:)) , cos(err(3,:)) );

% % the standard deviations for x,y,phi, computed by the filter:
% std_x = sqrt(P_sav(1,1,:));
% std_y = sqrt(P_sav(2,2,:));
% std_phi = sqrt(P_sav(3,3,:));

% find the estimated and corresponding true landmark positions:
% X_L_est = [X_k1k1(4:2:end)' ;X_k1k1(5:2:end)'] ;
weight_all = [particles.weight];
[w,idx] = max(weight_all);
X_L_est = particles(idx).Xl;
descriptors=particles(idx).dict;
X_L_true = X_L(:,descriptors);

%%% x-y plot of true and estimated trajectory
figure(1)
p1 = plot(X(1,:),X(2,:));
hold on
p2 = plot(Xr_sav(1,:),Xr_sav(2,:),'r');
pL = plot(X_L_true(1,:), X_L_true(2,:),'*');
pL_est = plot(X_L_est(1,:), X_L_est(2,:),'r*');
% plot the uncertainty ellipses for the landmarks
for i = 1:length(descriptors)
	plot_2D_ellipse(particles(idx).Xl(:,i), particles(idx).Pl(:,:,i), 0.99)
end
 
% add legend
if length(descriptors)>0
	legend([p1(1) p2(1) pL(1) pL_est(1)], 'True Robot Traj.', 'Est. Robot Traj', 'True Landmarks', 'Est. Landmarks')
else
	legend([p1(1) p2(1) ], 'True Robot Traj.', 'Est. Robot Traj')
end

% add axis labels
xlabel ('x (m)')
ylabel ('y (m)')

% plot of x-error and corresponding 3 sigma:
figure(2)
p1 = plot(err(1,:));
% hold on
% p2=plot(3*std_x(:), 'r');
% plot(-3*std_x(:), 'r');
% legend([p1(1) p2(1)], 'x error', '\pm 3 standard deviations')


% plot of y-error and corresponding 3 sigma:
figure(3)
p1 = plot(err(2,:));
% hold on
% p2=plot(3*std_y(:), 'r');
% plot(-3*std_y(:), 'r');
% legend([p1(1) p2(1)], 'y error', '\pm 3 standard deviations')


% plot of phi-error and corresponding 3 sigma:
figure(4)
p1 = plot(err(3,:));
% hold on
% p2=plot(3*std_phi(:), 'r');
% plot(-3*std_phi(:), 'r');
% legend([p1(1) p2(1)], '\phi error', '\pm 3 standard deviations')