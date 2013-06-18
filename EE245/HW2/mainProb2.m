clear all
close all

randn('seed',1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%parameters of the simulation:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% robot speed
v = 0.1; % m/sec
% robot rot. velocity
omega = 0.01; % rad/sec
% duration of simulation
timesteps = 1000;
% sampling interval for sensors
dt=1;
% std. deviation for odometry measurements:
sigma_v = 0.1; % m/sec
sigma_omega = 0.001; % rad/sec

% covariance matrix of GPS measurements
R_GPS = 1*eye(2); % m^2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%real world simulation:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
[X, v_m, omega_m, z] = rws(v, omega, timesteps, dt, sigma_v, sigma_omega, R_GPS);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the filter:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%%initial state estimate:
X_k1k1 = X(:,1);
P_k1k1 = zeros(3);

% save these for plotting:
X_sav(:,1) = X_k1k1;
P_sav(:,:,1) = P_k1k1;

for k=1:timesteps-1
 	
	% propagation
	[X_k1k ,P_k1k ] = ekf_propagate_2d(X_k1k1 ,P_k1k1, v_m(k),omega_m(k),dt,sigma_v,sigma_omega);

	% update
    [X_k1k1, P_k1k1] = ekf_update_2d(X_k1k,  P_k1k, z(:,k+1), R_GPS);
	
	% save these quantities
	X_sav(:,k+1) = X_k1k1(1:3);
	P_sav(:,:,k+1) = P_k1k1(1:3,1:3);
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

% compute the NEES for each timestep:  
for k = 1:timesteps
	nees(k) = err(:,k)'*inv(P_sav(:,:,k))*err(:,k);
end

%%% x-y plot of true and estimated trajectory
figure(1)
p1 = plot(X(1,:),X(2,:));
hold on
p2 = plot(X_sav(1,:),X_sav(2,:),'r'); 

% add legend
legend([p1(1) p2(1) ], 'True Robot Traj.', 'Est. Robot Traj')

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


 