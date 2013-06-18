function [theta_est, P_theta] = estimate_pose_GN(G_P, R_P, R, theta_init, N)

% iteration counter, to track progress
iter=0;

% the pose parameters:
G_P_R_i = theta_init(1:2);
phi_i = theta_init(3);

% a constant matrix we'll need....
J = [0 -1;
	1 0];

% maximum number of iterations
max_iter = 20;

while iter<max_iter

	% increase iteration counter
	iter=iter+1;

	% the estimated rotation matrix
	G_C_R_i = [cos(phi_i)  -sin(phi_i);
		         sin(phi_i)   cos(phi_i)];
	
    % initialization of the information matrix A and negative Jacobian vector b
    A = 0;
	b = 0;

	% loop through all points, and add the corresponding terms to A and b
    for j= 1:N
		% the measurement Jacobian
		H_j = G_C_R_i'*[-eye(2) -J*(G_P(:,j) - G_P_R_i)];
		
		% the estimated measurement (z_hat)
		z_hat_j = G_C_R_i'*(G_P(:,j) - G_P_R_i);
		
		% add the terms to build A and b
		A = A + H_j' *inv(R)* H_j;
		b = b + H_j'*inv(R)*(R_P(:,j) - z_hat_j);

	end
	
	% compute the correction to theta
	corr = inv(A)*b ;
	
	% update the estimates
	G_P_R_i = G_P_R_i + corr(1:2);
	phi_i = phi_i + corr(3);

	% stoppping condition
	if norm(corr)<1e-5
		break;
	end 
end

% the output estimate
theta_est = [G_P_R_i;phi_i];
% the covariance of the estimate:
P_theta = inv(A);
 
