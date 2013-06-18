function [X_k1k ,P_k1k ] = ekf_propagate_2d_SLAM(X_k1k1 ,P_k1k1, v_m,omega_m,dt,sigma_v,sigma_omega);

X_k1k = zeros(size(X_k1k1));
P_k1k = zeros(size(P_k1k1));

% assuming the upper 3 entries are X_r, the rest are X_Lk
% compute PhiR_k, equation 4.46
PhiR_k = [1.0   0.0   -v_m.*sin(X_k1k1(3))*dt;
         0.0   1.0   v_m.*cos(X_k1k1(3))*dt;
         0.0   0.0   1];

% compute GR_k, equation 4.47
GR_k = [ -cos(X_k1k1(3))*dt 0.0;
        -sin(X_k1k1(3))*dt 0.0;
        0.0           -1.0*dt];
    
% compute QR_k
QR_k = diag( [sigma_v sigma_omega] ).^2;
    
% compute X_k1k: using 4.42-4.44, 4.49
X_k1k(1) = X_k1k1(1) + v_m*cos(X_k1k1(3))*dt;
X_k1k(2) = X_k1k1(2) + v_m*sin(X_k1k1(3))*dt;
X_k1k(3) = X_k1k1(3) + omega_m*dt;
X_k1k(4:end) = X_k1k1(4:end);

% compute P_k1k, using 4.52
% Maybe compute Phi_k = [PhiR_k zeros; zeros eyes] is a better idea?
P_k1k(1:3,1:3) = PhiR_k*P_k1k1(1:3,1:3)*PhiR_k' + GR_k*QR_k*GR_k';
P_k1k(4:end,4:end) = P_k1k1(4:end,4:end);
P_k1k(1:3,4:end) = PhiR_k*P_k1k1(1:3,4:end);
P_k1k(4:end,1:3) = P_k1k1(4:end,1:3)*PhiR_k';