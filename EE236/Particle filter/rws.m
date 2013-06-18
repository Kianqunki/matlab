function [X, v_m, omega_m, z_GPS, z_th, z_r, X_L] = rws(v, omega, timesteps, dt, sigma_v, sigma_omega, R_GPS, sigma_theta,sigma_r)

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
X(:,1) = [0; -v/omega ; 0];

% compute the true robot poses:
for k=2:timesteps
    X(1,k) = X(1,k-1) + v_true(k-1)*dt*cos(X(3,k-1));
    X(2,k) = X(2,k-1) + v_true(k-1)*dt*sin(X(3,k-1));
    X(3,k) = X(3,k-1) + omega_true(k-1)*dt;
end

% create the GPS measurements:
z_GPS =zeros(2,timesteps);
sqrtmR = sqrtm(R_GPS);
for k = 1:timesteps
	z_GPS(:,k) = X(1:2,k) + sqrtmR * randn(2,1);
end

% create the range + bearing measurements

% choose the landmark position:
X_L = [v/omega/2 ; -v/omega/2 ];

z_th = zeros(1,timesteps);
z_r = zeros(1,timesteps);
for k = 1:timesteps
	% compute the landmark position in the robot frame (true one)
	C = [cos(X(3,k)) -sin(X(3,k)); 
	     sin(X(3,k))  cos(X(3,k))];
	R_P_L = C'*(X_L - X(1:2,k));
	
	% the range and bearing:
	[th,r] = cart2pol(R_P_L(1),R_P_L(2));
	z_th(k) = th + sigma_theta*randn;
	z_r(k) = r + sigma_r*randn;
end