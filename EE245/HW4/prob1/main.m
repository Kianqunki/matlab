close all
clear all  

load points_30_pct_out
%   load points_70_pct_out

% run estimation using all the measurements, as if no ouliers were present
[theta_est_plain, P_theta_plain] = estimate_displacement(R_p, G_p, R);

% run RANSAC
[theta_est_ransac, P_theta_ransac, inlier_indices_ransac] = estimate_displacement_ransac(R_p,G_p,R);

% run LMedS
[theta_est_lmeds, P_theta_lmeds, inlier_indices_lmeds] = estimate_displacement_lmeds(R_p,G_p,R);


% transform the measurements to the global frame using the estimates

% plain estimation
phi_est = theta_est_plain(3);
C_est = [cos(phi_est) -sin(phi_est);
	     sin(phi_est) cos(phi_est)];
G_p_est_plain = diag(theta_est_plain(1:2))*ones(2,length(R_p)) + C_est*R_p;

% RANSAC
phi_est = theta_est_ransac(3);
C_est = [cos(phi_est) -sin(phi_est);
	     sin(phi_est) cos(phi_est)];
G_p_est_ransac = diag(theta_est_ransac(1:2))*ones(2,length(inlier_indices_ransac)) + C_est*R_p(:,inlier_indices_ransac);
 
%LMedS
phi_est = theta_est_lmeds(3);
C_est = [cos(phi_est) -sin(phi_est);
	     sin(phi_est) cos(phi_est)];
G_p_est_lmeds = diag(theta_est_lmeds(1:2))*ones(2,length(inlier_indices_lmeds)) + C_est*R_p(:,inlier_indices_lmeds);
 

% plot the results.
figure
plot(G_p(1,:),G_p(2,:),'*')
hold on
plot(G_p_est_plain(1,:),G_p_est_plain(2,:),'c*')
plot(G_p_est_ransac(1,:),G_p_est_ransac(2,:),'g*')
axis equal
plot(G_p_est_lmeds(1,:),G_p_est_lmeds(2,:),'r*')
legend('Global features', 'Plain', 'RANSAC', 'LMedS')

% Saving the plot as .fig
hgsave('points_30_pct')
% Saving the plot as .jpg image
saveas(gcf,'points_30_pct.jpg');
 
% % Saving the plot as .fig
% hgsave('points_70_pct')
% % Saving the plot as .jpg image
% saveas(gcf,'points_70_pct.jpg');
