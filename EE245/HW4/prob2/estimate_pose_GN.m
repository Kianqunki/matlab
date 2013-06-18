function [theta_est, P_theta] = estimate_pose_GN(R1_P, R2_P, R1_cov, R2_cov, theta_init, N)

% iteration counter, to track progress
iter=0;

% the pose parameters:
R1_P_R2_i = theta_init(1:2);
phi_i = theta_init(3);

% a constant matrix we'll need....
J = [0 -1;
	1 0];


% infiinite loop, we will exit loop upon convergence 
while (iter < 10)

	% increase iteration counter
	iter=iter+1;

	% the estimated rotation matrix
	C_i = [cos(phi_i)  -sin(phi_i);
		         sin(phi_i)   cos(phi_i)];
	
    % initialization of the information matrix A and negative Jacobian vector b
    A = 0;
	b = 0;

	% loop through all points, and add the corresponding terms to A and b
    for j= 1:N
		% the measurement Jacobian
		H_j = [eye(2) J*C_i*R2_P(:,j)];
        H_n = [eye(2) -C_i];
        R = zeros(4);
        R(1:2,1:2) = R1_cov(:,:,j);
        R(3:4,3:4) = R2_cov(:,:,j);
        R = H_n*R*H_n';
		
		% the estimated measurement (z_hat)
		res = R1_P(:,j) - R1_P_R2_i - C_i*R2_P(:,j);
		
		% add the terms to build A and b
		A = A + H_j' *inv(R)* H_j;
		b = b + H_j'*inv(R)*res;

	end
	
	% compute the correction to theta
	corr = inv(A)*b ;
%     norm(corr)
	
	% update the estimates
	R1_P_R2_i = R1_P_R2_i + corr(1:2);
	phi_i = phi_i + corr(3);

	% stoppping condition
	if norm(corr)<1e-5
		break;
	end

end

% the output estimate
theta_est = [R1_P_R2_i;phi_i];
% the covariance of the estimate:
P_theta = inv(A);
 
