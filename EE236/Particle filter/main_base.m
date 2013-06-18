clear all
close all

% %set random seed
% seedRandom = 1;
% rand('state',seedRandom);
% randn('state',seedRandom);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%parameters of the simulation:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% robot speed
v = 0.1; % m/sec
% robot rot. velocity
omega = 0.01; % rad/sec
% duration of simulation
timesteps =2000;
% sampling interval for sensors
dt=1;
% std. deviation for odometry measurements:
sigma_v = 0.1; % m/sec
sigma_omega = 0.001; % rad/sec

% covariance matrix of GPS measurements
R_GPS = 0.2*eye(2);

% standard deviation of the range measurements
sigma_r = 0.1;% m
% standard deviation of the bearing measurements
sigma_theta = (2*pi/180);% rad


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filter parameters:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% number of particles
numParticle = 100;
% minimum allowed effective N
Nmin = 0.7*numParticle;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%real world simulation:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
[X, v_m, omega_m, z_GPS, z_theta, z_r, X_L] = rws(v, omega, timesteps, dt, sigma_v, sigma_omega, R_GPS, sigma_theta,sigma_r);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the filter:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%initialize the particles
% particles = repmat(struct('weight',1/numParticle,'Xr', [0;0;0],'Xl', [],'Pl',[]),numParticle, 1);
particles = repmat(struct('weight',1/numParticle,...
    'Xr', X(:,1)),...  % robot pose: position and orientation
    numParticle, 1);

% mean for plotting:
X_sav = zeros(3,timesteps);
X_sav(:,1) = X(:,1);

for k=1:timesteps-1	
	% plot timestep
	if mod(k,100)==0
		k
	end
	
    % propagation
    particles = pf_propagate_2d_SLAM(particles, v_m(k), omega_m(k), dt, sigma_v, sigma_omega);
    
    % update, with known data associtation
    particles = pf_update_2d_SLAM_GPS(particles, z_GPS(:,k+1), R_GPS);
%     particles = pf_update_2d_SLAM_range(particles, z_r(:,k+1), sigma_r, X_L);
%     particles = pf_update_2d_SLAM_bearing(particles, z_theta(:,k+1), sigma_theta, X_L);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Weight normalization and importance resampling
    particles = pf_resample( particles, Nmin );
     	
% 	% save max weight particles for plotting
%     weight_all = [particles.weight];
%     [w,idx] = max(weight_all);
%     X_sav(:,k+1) = particles(idx).Xr(:);
    
	% save mean for plotting
	X_sav(:,k+1) = find_pf_mean(particles);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotting:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%%% x-y plot of true and estimated trajectory
figure(1)
p1 = plot(X(1,:),X(2,:));
hold on
p2 = plot(X_sav(1,:),X_sav(2,:),'r');
poses = [particles.Xr];
plot(poses(1,:),poses(2,:),'.')
plot(X_L(1),X_L(2),'*')
% add legend
legend([p1(1) p2(1)], 'True', 'Estimate')
% add axis labels
xlabel ('x (m)')
ylabel ('y (m)') 