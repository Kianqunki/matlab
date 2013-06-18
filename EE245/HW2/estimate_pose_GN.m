function [theta_est, P_theta] = estimate_pose_GN(G_p, R_p, R, theta_init, N)

% starting Gauss-Newton
max_iteration = 20;
convergence_norm = 1e-8;

theta_est = theta_init;
P_theta = zeros(3);
J = [0 -1; 1 0];

for l=1:max_iteration
    A = zeros(3,3);
    b = zeros(3,1);
    
    pos_hat = theta_est(1:2,1);
    phi_hat = theta_est(3);
    G_C_R = [cos(phi_hat) -sin(phi_hat); sin(phi_hat) cos(phi_hat)];
    
    % for all measurements
    for i=1:N
        z_hat = G_C_R'*(G_p(:,i) - pos_hat);
        % residual
        r = R_p(:,i) - z_hat;
        
        % measurement Jacobian
        H_i = zeros(2,3);
        H_i(:,1:2) = -G_C_R';
        H_i(:,3) = -G_C_R'*J*(G_p(:,i) - pos_hat);
        
        % compute Jacobian of the cost function
        b = b - H_i'*inv(R)*r;
        % compute Hessian of the cost function
        A = A + H_i'*inv(R)*H_i;
    end
    
    % compute the correction
    delta_theta = -A\b;
    theta_est = theta_est + delta_theta;
    P_theta = inv(A);
    
    norm(delta_theta)
    
    if ( norm(delta_theta) < convergence_norm )
        break
    end
    
end

end