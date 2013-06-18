function [X,v_m,omega_m,z,X_L] = rws(v, omega, timesteps, dt, sigma_v, sigma_omega, sigma_th, min_range, max_range, L_N, r_p_s1, r_p_s2, r_p_s3)

noise = 1;

% the true velocities for the duration of the simulation:
v_true = v*ones(1,timesteps);
omega_true = omega*ones(1,timesteps);

%velocity measurements
v_m = v_true + noise*sigma_v*randn(1,timesteps);
omega_m = omega_true + noise*sigma_omega*randn(1,timesteps);

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
sqrtmR = sigma_th*eye(2); 
for k = 1:timesteps
	% the distances of all landmarks to the robot
	d = sqrt((X_L(1,:)-X(1,k)).^2 + (X_L(2,:)-X(2,k)).^2 );
	% the indices to the landmarks we can measure
	ind = find(d<max_range & d>min_range);
	
	% create the measurements
	
	% first find the positions of the sensors in G:
	C = [cos(X(3,k)) -sin(X(3,k));
		 sin(X(3,k))  cos(X(3,k))];
	G_p_s1 = X(1:2,k) + C*r_p_s1;
	G_p_s2 = X(1:2,k) + C*r_p_s2; 
	
	% get the vectors from the sensors to the landmarks we see, in the S
	% frames
	v1 = C' * (X_L(:,ind) - diag(G_p_s1)*ones(2,length(ind)));
	v2 = C' * (X_L(:,ind) - diag(G_p_s2)*ones(2,length(ind))); 
	
	% get the bearing angles to the landmarks:
	th1 = cart2pol(v1(1,:),v1(2,:));
	th2 = cart2pol(v2(1,:),v2(2,:)); 
	
	% only keep measurements to landmarks that have sufficient "baseline"
	good_ind = find(abs(th1-th2)>4*pi/180);
	th1 = th1(good_ind);
	th2 = th2(good_ind);
	ind = ind(good_ind);
	
	if length(ind)==0
		z(k).measurements = [];
		z(k).descriptors = [];
		continue;
	end
	
	% get the measurements:
	z_k = [th1;th2] + noise*sqrtmR*randn(2,length(ind));
	z(k).measurements = z_k;
	z(k).descriptors = ind;
end