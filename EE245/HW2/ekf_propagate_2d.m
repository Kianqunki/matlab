function [X_k1k ,P_k1k ] = ekf_propagate_2d(X_k1k1 ,P_k1k1, v_m, omega_m, dt, sigma_v, sigma_omega);

X_k1k = zeros(size(X_k1k1));
P_k1k = zeros(size(P_k1k1));

% compute PhiR_k
PhiR_k = [1.0   0.0   -v_m.*sin(X_k1k1(3))*dt;
         0.0   1.0   v_m.*cos(X_k1k1(3))*dt;
         0.0   0.0   1];

% compute GR_k
GR_k = [ -cos(X_k1k1(3))*dt 0.0;
        -sin(X_k1k1(3))*dt 0.0;
        0.0           -1.0*dt];
    
% compute QR_k
QR_k = diag( [sigma_v sigma_omega] ).^2;
    
% compute X_k1k
X_k1k(1) = X_k1k1(1) + v_m*cos(X_k1k1(3))*dt;
X_k1k(2) = X_k1k1(2) + v_m*sin(X_k1k1(3))*dt;
X_k1k(3) = X_k1k1(3) + omega_m*dt;

% compute P_k1k
P_k1k(1:3,1:3) = PhiR_k*P_k1k1(1:3,1:3)*PhiR_k' + GR_k*QR_k*GR_k';