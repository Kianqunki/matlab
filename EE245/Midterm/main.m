clear all
close all

randn('seed',133)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%parameters of the simulation:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% robot speed
v = 0.2; % m/sec
% robot rot. velocity
omega = 0.01; % rad/sec
% duration of simulation
timesteps = 2000;
% sampling interval for sensors
dt=1;
% std. deviation for odometry measurements:
sigma_v = 0.02; % m/sec
sigma_omega = 0.002; % rad/sec

% the positions of the sensors on the robot
r_p_s1 =[1;0];
r_p_s2 =[1;1]; 

% std of the bearing measurements
sigma_th = 0.5*pi/180; % 0.5 degrees, expressed in rad

% sensors range 
max_range = 6; %m
min_range = 2; %m

% we'll create landmarms on an L_N x L_N grid (not all landmarks will be within
% the robot's sensing range though)
L_N = 9; % we create a pair of landmarks in each location

% set this to true to use the robot-to-landmark measurements, or to false to use
% the odometry only
use_updates = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%real world simulation:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
[X, v_m, omega_m, z, X_L] = rws_2d_SLAM(v, omega, timesteps, dt, sigma_v, sigma_omega,  sigma_th, min_range, max_range, L_N,r_p_s1,r_p_s2);
params.X = X;
params.X_L = X_L;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the filter:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%%initial state estimate:
X_k1k1 = X(:,1);
P_k1k1 = zeros(3);

% we also maintain a vector with the descriptors of the landmarks that we
% have added to the state vector:
descriptors = [];  % initially we don't have any landmarks

% save these for plotting:
X_sav(:,1) = X_k1k1;
P_sav(:,:,1) = P_k1k1;

for k=1:timesteps-1
 	
	% propagation
	[X_k1k ,P_k1k ] = ekf_propagate_2d_SLAM(X_k1k1 ,P_k1k1, v_m(k),omega_m(k),dt,sigma_v,sigma_omega);

	% update
	if use_updates == true
		[X_k1k1, P_k1k1, descriptors] = ekf_update_2d_SLAM(X_k1k,  P_k1k, z(k+1), sigma_th, r_p_s1, r_p_s2, descriptors, params);
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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotting:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% compute the estimation error for the robot:
err = X - X_sav;
% make sure the orientation error is in the first cycle:
err(3,:) = atan2( sin(err(3,:)) , cos(err(3,:)) );

% the standard deviations for x,y,phi, computed by the filter:
std_x = sqrt(P_sav(1,1,:));
std_y = sqrt(P_sav(2,2,:));
std_phi = sqrt(P_sav(3,3,:));

% find the estimated and corresponding true landmark positions:
X_L_est = [X_k1k1(4:2:end)' ;X_k1k1(5:2:end)'] ;
X_L_true = X_L(:,descriptors);

% compute the NEES for each timestep:  
for k = 1:timesteps
	nees(k) = err(:,k)'*inv(P_sav(:,:,k))*err(:,k);
end


%%% x-y plot of true and estimated trajectory
figure(1)
p1 = plot(X(1,:),X(2,:));
hold on
p2 = plot(X_sav(1,:),X_sav(2,:),'r');
pL = plot(X_L_true(1,:), X_L_true(2,:),'*');
pL_est = plot(X_L_est(1,:), X_L_est(2,:),'r*');
% plot the uncertainty ellipses for the landmarks
for i = 1:size(X_L_est,2)
	plot_2D_ellipse(X_L_est(:,i), P_k1k1(3+2*i-1:3+2*i,3+2*i-1:3+2*i), 0.99)
end
 

% add legend
if length(descriptors)>0
	legend([p1(1) p2(1) pL(1) pL_est(1)], 'True Robot Traj.', 'Est. Robot Traj', 'True Land.', 'Est. Land.')
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
p2=plot(3*std_x(:), 'r');
plot(-3*std_x(:), 'r');
legend([p1(1) p2(1)], 'x error', '\pm 3 standard deviations')


% plot of y-error and corresponding 3 sigma:
figure(3)
p1 = plot(err(2,:));
hold on
p2=plot(3*std_y(:), 'r');
plot(-3*std_y(:), 'r');
legend([p1(1) p2(1)], 'y error', '\pm 3 standard deviations')


% plot of phi-error and corresponding 3 sigma:
figure(4)
p1 = plot(err(3,:));
hold on
p2=plot(3*std_phi(:), 'r');
plot(-3*std_phi(:), 'r');
legend([p1(1) p2(1)], '\phi error', '\pm 3 standard deviations')

% plot NEES
figure(5)
hold on
plot(nees)
title('NEES of robot pose')


 