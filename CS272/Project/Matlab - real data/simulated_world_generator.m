function [X,v_m,omega_m,z,X_L] = rws(v, omega, timesteps, dt, sigma_v, sigma_omega, R, min_range, max_range, L_N)
 
% the true velocities for the duration of the simulation:
v_true = v*ones(1,timesteps);
omega_true = omega*ones(1,timesteps);

%velocity measurements
v_m = v_true + sigma_v*randn(1,timesteps);
omega_m = omega_true + sigma_omega*randn(1,timesteps);

% compute the true robot pose

% initialize the matrix (this speeds up the code, by allocating memory 
% just once, and not every time a new column is added to the matrix. 
% It can make a huge difference in runtime when the number of 
% timesteps is large):
X = zeros(3,timesteps);

% compute the true robot poses:
for k=2:timesteps
    X(1,k) = X(1,k-1) + v_true(k-1)*dt*cos(X(3,k-1));
    X(2,k) = X(2,k-1) + v_true(k-1)*dt*sin(X(3,k-1));
    X(3,k) = X(3,k-1) + omega_true(k-1)*dt;
end

% create the landmark positions:
% find the radius of the trajectory
radius = v/omega;
% get the landmark placement
a= linspace(-radius-max_range, radius+max_range, L_N );
[x_l, y_l] = meshgrid(a);
X_L = [x_l(:)'; radius+y_l(:)']; 


% now get the measurements
sqrtmR = sqrtm(R); 
for k = 1:timesteps
	% the distances of all landmarks to the robot
	d = sqrt((X_L(1,:)-X(1,k)).^2 + (X_L(2,:)-X(2,k)).^2 );
	% the indices to the landmarks we can measure
	ind = find(d<max_range & d>min_range);
	
	% create the measurements:
	C = [cos(X(3,k)) -sin(X(3,k));
		 sin(X(3,k))  cos(X(3,k))];
	z_k = C' * (X_L(:,ind) - diag(X(1:2,k))*ones(2,length(ind))) + sqrtmR*randn(2,length(ind));
	z(k).measurements = z_k;
	z(k).descriptors = ind;
end