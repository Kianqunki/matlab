function [theta_est, P_theta, inlier_indices] = estimate_displacement_lmeds(R_p, G_p, R)

worst_inlier_fraction = 0.3;
% The maximum number of trials (sampling)
max_trial_no = ceil(log(1-0.99)/log(1-worst_inlier_fraction^2));
% Total number of feature points
point_no = size(R_p,2);
% Thresholde for Mahalanobis distance gating test
thres = chi2inv(0.9,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get initial guess:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
best_score = 1e8;
best_theta = zeros(3,1);

for i=1:max_trial_no
    % generate two random index between 1 and point_no
    sample_index = ceil(1 + (point_no-1)*rand(2,1));
    if ( sample_index(1) == sample_index(2) )
        sample_index(2) = sample_index(1) - 1;
    end
    
    % the difference of the vectors in the robot frame:
    R_dp= R_p(:,sample_index(1)) - R_p(:,sample_index(2));
    % the difference of the vectors in the global frame:
    G_dp= G_p(:,sample_index(1)) - G_p(:,sample_index(2));
    
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
    G_p_R_init = G_p(:,sample_index(1)) - C_init*R_p(:,sample_index(1));
    
    % the initial estimate for G-N minimization:
    theta_init = [G_p_R_init ;phi_init];

    % Compute Mahalanobid distance
    maha_dist = zeros(point_no,1);
    for j=1:point_no
        diff = R_p(:,j) + C_init'*(G_p_R_init - G_p(:,j));
        maha_dist(j) = 1/R*diff'*diff;
    end
    
    % Find points that pass the Mahalanobis distance gating test
    idx = find(maha_dist < thres);
    score = median(maha_dist);
    
    if ( score < best_score )
        best_score = score;
        best_theta = theta_init;
        inlier_indices = idx;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run Gauss-Newton:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[theta_est, P_theta] = estimate_pose_GN(G_p(:,inlier_indices), R_p(:,inlier_indices), R, best_theta, length(inlier_indices));