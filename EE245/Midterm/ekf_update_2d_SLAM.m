function [X_k1k1, P_k1k1, descriptors] = ekf_update_2d_SLAM(X_k1k,  P_k1k, z, sigma_th, r_p_s1, r_p_s2, descriptors, params)

J = [0 -1; 1 0];
X_k1k1 = X_k1k;
P_k1k1 = P_k1k;
R = sigma_th^2*eye(2);

for i=1:length(z.descriptors)
    idx = find( descriptors == z.descriptors(i),1);
    
    % if the landmark is already discovered
    if ( length(idx) == 1)
        C = [cos(X_k1k1(3)) -sin(X_k1k1(3)); 
            sin(X_k1k1(3)) cos(X_k1k1(3))];
        X_R = X_k1k1(1:2);
        X_L = X_k1k1(2*idx+2:2*idx+3);   % estimated location of landmark
        
        z_k1 = z.measurements(:,i);
        
        % find bearing measurement from sensor 1
        q1 = C'*(X_L - X_R) - r_p_s1;
        t1 = q1(2)/q1(1);
        z1j = atan2(q1(2),q1(1));
        % find bearing measurement from sensor 2
        q2 = C'*(X_L - X_R) - r_p_s2;
        t2 = q2(2)/q2(1);
        z2j = atan2(q2(2),q2(1));
        
        z_hat_k1 = [z1j; z2j];
        r = z_k1 - z_hat_k1;
%         X_k1k1
%         keyboard;
        
        % compute Jacobian
        HR_k1 = [-C' -C'*J*(X_L - X_R)];
        HL_k1 = C';
        Delta1 = 1/(1+t1^2)*[-q1(2)/q1(1)^2 1/q1(1)];
        Delta2 = 1/(1+t2^2)*[-q2(2)/q2(1)^2 1/q2(1)];
        
        H_k1 = zeros(2, length(X_k1k1));
        H_k1(1,1:3) = Delta1*HR_k1;
        H_k1(2,1:3) = Delta2*HR_k1;
        H_k1(1,2*idx+2:2*idx+3) = Delta1*HL_k1;
        H_k1(2,2*idx+2:2*idx+3) = Delta2*HL_k1;
        
        % EKF update equations
        S_k = H_k1*P_k1k1*H_k1' + R;
        %K_k1 = P_k1k1*H_k1'*inv(S_k);
        K_k1 = P_k1k1*H_k1'/(S_k);
        X_k1k1 = X_k1k1 + K_k1*r;
        P_k1k1 = (eye(length(P_k1k1)) - K_k1*H_k1)*P_k1k1*(eye(length(P_k1k1)) - K_k1*H_k1)' + K_k1*R*K_k1';
    else % if the landmark is not discovered
        C = [cos(X_k1k1(3)) -sin(X_k1k1(3));
            sin(X_k1k1(3)) cos(X_k1k1(3))];
        X_R = X_k1k1(1:2);
        z_k1 = z.measurements(:,i);
        a1 = z_k1(1);
        a2 = z_k1(2);
        
        % initialize the landmark
        phi = (r_p_s1 + 1/sin(a1 - a2)*[sin(a2) -cos(a2)]*(r_p_s1-r_p_s2)*[cos(a1); sin(a1)]);
        X_L = C*phi + X_R;
%         keyboard;
        Theta1 = [eye(2) J*C*phi];
        G = zeros(2);
        alpha = a1 - a2;
        diff_p = (r_p_s1-r_p_s2);
        G(:,1) = 1/sin(alpha)*[sin(a2) -cos(a2)]*diff_p*[sin(a1); -cos(a1)] + cos(alpha)/sin(alpha)^2*[sin(a2) -cos(a2)]*diff_p*[cos(a1); sin(a1)];
        G(:,2) = -cos(alpha)/sin(alpha)^2*[sin(a2) -cos(a2)]*diff_p*[cos(a1); sin(a1)] + 1/sin(alpha)*[-cos(a2) -sin(a2)]*diff_p*[cos(a1); sin(a1)];
        Theta2 = C*G;
        
        % augment the state vector and the covariance matrix
        entryNum = size(X_k1k1,1);
        X_k1k1 = [X_k1k1; X_L];
        Prr_k1 = P_k1k1(1:3,1:3);
        P_k1k1_new = zeros(entryNum+2);
        P_k1k1_new(1:entryNum,1:entryNum) = P_k1k1;
        
        % P_LnR and P_LnLi
        P_k1k1_new(entryNum+1:entryNum+2,1:entryNum) = Theta1*P_k1k1(1:3,1:entryNum);
        % P_RLn and P_LiLn
        P_k1k1_new(1:entryNum,entryNum+1:entryNum+2)= P_k1k1(1:entryNum,1:3)*Theta1';
        P_k1k1_new(entryNum+1:entryNum+2,entryNum+1:entryNum+2) = Theta1*Prr_k1*Theta1'+Theta2*R*Theta2';
        
        %update the output
        P_k1k1 = P_k1k1_new;
        descriptors = [descriptors z.descriptors(i)];
    end
end

end