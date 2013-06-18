function [X_k1k1, P_k1k1] = ekf_update_2d(X_k1k, P_k1k, z, R_GPS)

% Compute C^T(phi_R)
J = [0 -1; 1 0];
X_k1k1 = X_k1k;
P_k1k1 = P_k1k;

% the actual measurement
z_k1 = z;
% the predicted measurement
z_hat_k1 = X_k1k1(1:2);
% measurement residual 
r = z_k1 - z_hat_k1;
% the measurement Jacobian
H_k1 = zeros(2,3);
H_k1(1:2,1:2) = eye(2);

% the measurement covariance
S_k = H_k1*P_k1k1*H_k1' + R_GPS;
% the Kalman gain
%K_k1 = P_k1k1*H_k1'*inv(S_k);
K_k1 = P_k1k1*H_k1'/(S_k);
% the updated state
X_k1k1 = X_k1k1 + K_k1*r;
% the updated covariance matrix
P_k1k1 = P_k1k1 - K_k1*H_k1*P_k1k1;