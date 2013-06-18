clear all
close all

randn('seed',1)
DATA_KNOWN = 0;

% setting up the simulation

% Parameters for process model and measurement model
% robot's velocity v
vel = 0.3; % m/sec
% robot's rotational velocity omega
omega = 0.015; % rad/sec
% time length of simulation, aka number of timesteps
timesteps = 2000;
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
%[X, v_m, omega_m, z, X_L] = simulated_world_generator(vel, omega, timesteps, dt, sigma_vel, sigma_omega, R, min_range, max_range, L_N);
% real data from http://www-personal.acfr.usyd.edu.au/nebot/car_park.htm
[X, v_m, omega_m, dt, z, X_L] = real_world_generator();

% NOTE: tried to generate the real world data to fit with defined format by
% simulator. However, the measurement model and process model are much different and need to be
% revised.
% aka reprogramming the whole propagation and update function

% Currently: Process model:
% Steering propagation is different
% ...
% Measurement model:
% laser is not at {R} robot coordinate frame
% Measurement as distance and bearing, not relative position (convertible
% though)

timesteps = length(v_m);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the filter: 

%%initial state estimate:
X_k1k1 = [X(1,1);X(2,1);0];
P_k1k1 = zeros(3);

% we also maintain a vector with the descriptors of the landmarks that we
% have added to the state vector:
descriptors = [];  % initially we don't have any landmarks
dict=[];

% save the robot positions for plotting trajectory:
Xr_sav = zeros(length(X_k1k1),timesteps);
Pr_sav(:,:,1) = P_k1k1;

for k=1:timesteps-1
 	
	% propagation
	%[X_k1k ,P_k1k ] = ekf_propagate_2d_SLAM(X_k1k1 ,P_k1k1, v_m(k),omega_m(k),dt,sigma_vel,sigma_omega);
    [X_k1k ,P_k1k ] = ekf_propagate_2d_SLAM(X_k1k1 ,P_k1k1, v_m(k),omega_m(k),dt(k),sigma_vel,sigma_omega);

	% update
    if (DATA_KNOWN)
        [X_k1k1, P_k1k1, descriptors] = ekf_update_2d_SLAM(X_k1k,  P_k1k, z(k+1), R, descriptors);
    else
        if (~isempty(z(k+1).measurements))
        %if (false)
            [X_k1k1, P_k1k1, descriptors, dict] = ekf_update_associate_2d_SLAM(X_k1k,  P_k1k, z(k+1), R, descriptors, dict);
        else
            X_k1k1 = X_k1k;
            P_k1k1 = P_k1k;
        end
    end
		
	% save these quantities for plotting
	Xr_sav(:,k+1) = X_k1k1(1:3);
	Pr_sav(:,:,k+1) = P_k1k1(1:3,1:3);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotting:

% % compute the estimation error for the robot:
% err = X - Xr_sav;
% % make sure the orientation error is in the first cycle:
% err(3,:) = atan2( sin(err(3,:)) , cos(err(3,:)) );
% 
% % the standard deviations for robot pose:
% std_x = sqrt(Pr_sav(1,1,:));
% std_y = sqrt(Pr_sav(2,2,:));
% std_phi = sqrt(Pr_sav(3,3,:));

% find the estimated and corresponding true landmark positions:
X_L_est = [X_k1k1(4:2:end)' ;X_k1k1(5:2:end)'] ;
if (DATA_KNOWN)
    X_L_true = X_L(:,descriptors);
else
    X_L_true = X_L(:,dict);
end


%%% x-y plot of true and estimated trajectory
figure(1)
p1 = plot(X(1,:),X(2,:));
hold on
p2 = plot(Xr_sav(1,:),Xr_sav(2,:),'r');
pL = plot(X_L_true(1,:), X_L_true(2,:),'*');
pL_est = plot(X_L_est(1,:), X_L_est(2,:),'r*');
% plot the uncertainty ellipses for the landmarks
for i = 1:length(descriptors)
	plot_2D_ellipse(X_L_est(:,i), P_k1k1(3+2*i-1:3+2*i,3+2*i-1:3+2*i), 0.99)
%  	plot_2D_ellipse_t(X_L_est(:,i), P_k1k1(3+2*i-1:3+2*i,3+2*i-1:3+2*i), 0.99)
end

% % add legend
% if length(descriptors)>0
% 	legend([p1(1) p2(1) pL(1) pL_est(1)], 'True Robot Traj.', 'Est. Robot Traj', 'True Landmarks', 'Est. Landmarks')
% else
% 	legend([p1(1) p2(1) ], 'True Robot Traj.', 'Est. Robot Traj')
% end

% add axis labels
xlabel ('x (m)')
ylabel ('y (m)')

% % plot of x-error and corresponding 3 sigma:
% figure(2)
% p1 = plot(err(1,:));
% hold on
% p2=plot(3*std_x(:), 'r');
% plot(-3*std_x(:), 'r');
% legend([p1(1) p2(1)], 'x error', '\pm 3 standard deviations')
% 
% 
% % plot of y-error and corresponding 3 sigma:
% figure(3)
% p1 = plot(err(2,:));
% hold on
% p2=plot(3*std_y(:), 'r');
% plot(-3*std_y(:), 'r');
% legend([p1(1) p2(1)], 'y error', '\pm 3 standard deviations')
% 
% 
% % plot of phi-error and corresponding 3 sigma:
% figure(4)
% p1 = plot(err(3,:));
% hold on
% p2=plot(3*std_phi(:), 'r');
% plot(-3*std_phi(:), 'r');
% legend([p1(1) p2(1)], '\phi error', '\pm 3 standard deviations')