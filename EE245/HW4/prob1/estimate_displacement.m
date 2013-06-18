function [theta_est, P_theta] = estimate_displacement(R_p, G_p, R)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% get initial guess:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% the difference of the vectors in the robot frame:
R_dp= R_p(:,1) - R_p(:,end);
% the difference of the vectors in the global frame:
G_dp= G_p(:,1) - G_p(:,end);

%build the A matrix of the linear system:
A =[R_dp(1)  -R_dp(2);
	R_dp(2)    R_dp(1)];

% the vector v
v = inv(A)*G_dp;

% the orientation estimate:
phi_init = atan2(v(2), v(1));

% the rotation matrix corresponding to phi_init
C_init = [cos(phi_init) -sin(phi_init);
	      sin(phi_init) cos(phi_init)];

% the position estimate 
G_p_R_init = G_p(:,1) - C_init*R_p(:,1);

% the initial estimate for G-N minimization:
theta_init = [G_p_R_init ;phi_init];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% run Gauss-Newton:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[theta_est, P_theta] = estimate_pose_GN(G_p, R_p, R, theta_init,size(R_p,2)); 