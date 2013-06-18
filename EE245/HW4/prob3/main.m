clear all
close all

randn('seed',19)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%parameters of the simulation:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% robot speed
v = 0.3; % m/sec
% robot rot. velocity
omega = 0.015; % rad/sec
% duration of simulation
timesteps = 1400;
dt =1;

% std. deviation for odometry measurements:
sigma_dx = 0.04; % m/sec
sigma_dy = 0.03; % m/sec
sigma_dphi = 0.7*pi/180; % rad/sec
Q = diag([sigma_dx^2;sigma_dy^2;sigma_dphi^2]);

% covariance matrix of landmark measurements
R = 1*eye(2); % m^2

% prior covariance (covariance of estimate for first pose)
P_prior = 1e-3*eye(3);

% range of laser scanner
max_range = 6; %m
min_range = 0.4; %m

% we'll create landmarms on an L_N x L_N grid (not all landmarks will be within
% the robot's sensing range though)
L_N = 9; % we create a pair of landmarks in each location

% set this to true to use the robot-to-landmark measurements, or to false to use
% the odometry only
use_updates = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%real world simulation:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
[X, dp_meas, z, X_L] = rws_2d_SLAM(v, omega, timesteps, dt, Q, R, min_range, max_range, L_N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the filter:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%%initial state estimate:
X_1_prior = X(:,1)+ sqrtm(P_prior)*randn(3,1);
X_k1k1 = X_1_prior;
P_k1k1 = P_prior;

% we also maintain a vector with the descriptors of the landmarks that we
% have added to the state vector:
descriptors = [];  % initially we don't have any landmarks

% save these for plotting:
X_sav(:,1) = X_k1k1;
P_sav(:,:,1) = P_k1k1;

for k=1:timesteps-1
 	
	% propagation
	[X_k1k ,P_k1k ] = ekf_propagate_2d_SLAM(X_k1k1 ,P_k1k1, dp_meas(:,k),dt,Q);

	% update
	if use_updates == true
		[X_k1k1, P_k1k1, descriptors] = ekf_update_2d_SLAM(X_k1k,  P_k1k, z(k+1), R, descriptors);
	else
		X_k1k1 = X_k1k;
		P_k1k1 = P_k1k;
	end
	
	% save these quantities, we'll need them for plotting
	X_sav(:,k+1) = X_k1k1(1:3);
	P_sav(:,:,k+1) = P_k1k1(1:3,1:3);
	
	%printout so we know things are moving on
	if mod(k,200)==0
		k
	end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get the initial guess for the state used in MAP:
N_landmarks = (length(X_k1k1)-3)/2;
X_init =zeros(3*timesteps+2*N_landmarks,1);
for i = 1:timesteps
	X_init(3*i-2:3*i) = X_sav(:,i);
end
X_init(3*timesteps+1:end) = X_k1k1(4:end);

% run the MAP estimator
X_est_map = map_estimate(X_init, z, dp_meas, descriptors, X_1_prior, P_prior, timesteps, Q, R);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotting:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% the MAP estimate for the trajectory
X_sav_map = reshape(X_est_map(1:3*timesteps),3,timesteps);

% compute the estimation error for the robot:
err = X - X_sav;
err_map = X - X_sav_map;
% make sure the orientation error is in the first cycle:
err(3,:) = atan2( sin(err(3,:)) , cos(err(3,:)) );
err_map(3,:) = atan2( sin(err_map(3,:)) , cos(err_map(3,:)) );

% the standard deviations for x,y,phi, computed by the EKF:
std_x = sqrt(P_sav(1,1,:));
std_y = sqrt(P_sav(2,2,:));
std_phi = sqrt(P_sav(3,3,:));

% find the estimated and corresponding true landmark positions:
X_L_est = [X_k1k1(4:2:end)' ;X_k1k1(5:2:end)'] ;
X_L_est_map = reshape(X_est_map(3*timesteps+1:end),2,size(X_L_est,2));
X_L_true = X_L(:,descriptors);
 
%%% x-y plot of true and estimated trajectory
figure(1)
p1 = plot(X(1,:),X(2,:));
hold on
p2 = plot(X_sav(1,:),X_sav(2,:),'r');
p2_map = plot(X_sav_map(1,:),X_sav_map(2,:),'c');
pL = plot(X_L_true(1,:), X_L_true(2,:),'*');
pL_est = plot(X_L_est(1,:), X_L_est(2,:),'r*');
pL_est_map = plot(X_L_est_map(1,:), X_L_est_map(2,:),'g*');
% plot the uncertainty ellipses for the landmarks
for i = 1:size(X_L_est,2)
	plot_2D_ellipse(X_L_est(:,i), P_k1k1(3+2*i-1:3+2*i,3+2*i-1:3+2*i), 0.99)
end
 

% add legend
if length(descriptors)>0
	legend([p1(1) p2(1) p2_map(1) pL(1) pL_est(1) pL_est_map(1)], 'True Robot Traj.', 'EKF Robot Traj', 'MAP Robot Traj' ,'True Land.', 'EKF Land.', 'MAP Land.')
else
	legend([p1(1) p2(1) ], 'True Robot Traj.', 'Est. Robot Traj')
end

% add axis labels
xlabel ('x (m)')
ylabel ('y (m)')



% plot of x-error and corresponding 3 sigma:
figure(2)
p1 = plot(err(1,:));
hold on
p1_map = plot(err_map(1,:),'c');
p2=plot(3*std_x(:), 'r');
plot(-3*std_x(:), 'r');
legend([p1(1)  p1_map(1) p2(1)], 'x error', 'x error MAP','\pm 3 standard deviations')


% plot of y-error and corresponding 3 sigma:
figure(3)
p1 = plot(err(2,:));
hold on
p1_map = plot(err_map(2,:),'c');
p2=plot(3*std_y(:), 'r');
plot(-3*std_y(:), 'r');
legend([p1(1) p1_map(1) p2(1)], 'y error', 'y error MAP', '\pm 3 standard deviations')


% plot of phi-error and corresponding 3 sigma:
figure(4)
p1 = plot(err(3,:));
hold on
p1_map = plot(err_map(3,:),'c');
p2=plot(3*std_phi(:), 'r');
plot(-3*std_phi(:), 'r');
legend([p1(1) p1_map(1) p2(1)], '\phi error', '\phi error MAP','\pm 3 standard deviations')
 
 