clear all
close all

load points2D.mat

% measurement covariance matrix
R = 0.04^2*eye(2);

% find the initial guess
R_p_1 = R_p(:,1);
R_p_N = R_p(:,N);
G_p_1 = G_p(:,1);
G_p_N = G_p(:,N);

% the vector b
b = G_p_1 - G_p_N;
% the matrix A
A = zeros(2);
A(1,1) = R_p_1(1) - R_p_N(1);
A(1,2) = -(R_p_1(2) - R_p_N(2));
A(2,1) = (R_p_1(2) - R_p_N(2));
A(2,2) = A(1,1);

v = A\b;
init_phi = atan2(v(2), v(1));

rot_matrix = [cos(init_phi) -sin(init_phi); sin(init_phi) cos(init_phi)];
init_pos = G_p_1 - rot_matrix*R_p_1;
theta_init = [init_pos; init_phi];

% Question ii
[theta_est, P_theta] = estimate_pose_GN(G_p, R_p, R, theta_init, N);

theta_est

P_theta

% Question iii
[theta_est, P_theta] = estimate_pose_GN(G_p, R_p, R, theta_init, 1);