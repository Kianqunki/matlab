function [X,v_m,omega_m,z] = rws(v, omega, timesteps, dt, sigma_v, sigma_omega, R_GPS)
 
% the true velocities for the duration of the simulation:
v_true = v*ones(1,timesteps);
omega_true = omega*ones(1,timesteps);

%velocity measurements
v_m = v_true + sigma_v*randn(1,timesteps);
omega_m = omega_true + sigma_omega*randn(1,timesteps);

% compute the true robot pose

X = zeros(3, timesteps);
z = zeros(2, timesteps);

% compute the true robot poses:
for k=2:timesteps
    X(1,k) = X(1,k-1) + v_true(k-1)*dt*cos(X(3,k-1));
    X(2,k) = X(2,k-1) + v_true(k-1)*dt*sin(X(3,k-1));
    X(3,k) = X(3,k-1) + omega_true(k-1)*dt;
end

% now get the measurements
sqrtmR = sqrtm(R_GPS); 
for k = 1:timesteps	
	% create the measurements:
	z(:,k) = X(1:2,k) + sqrtmR*randn(2,1);
end