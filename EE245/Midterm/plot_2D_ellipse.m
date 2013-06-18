function plot_2D_ellipse(m, C, p)
% m : the 2x1 mean
% C : the 2x2 covariance matrix
% p : the probability level, in [0,1]

% create N points at regular angle intervals:
N = 400;

% the angles
theta_i = linspace(0,2*pi, N);

% the points
p_i = [cos(theta_i);sin(theta_i)];

% the parameter defining the size of the ellipse:
gamma = chi2inv(p,2);

% we need to compute the norm of the vector at each of these orientations.
for i = 1:length(theta_i)
	rho_i(i) = sqrt(gamma/(p_i(:,i)'*inv(C)*p_i(:,i)));
end
	
% the points on the ellipse:
for i =1:length(theta_i)
	p_ell(:,i) = m + rho_i(i)*p_i(:,i);
end

% plot the points
plot(p_ell(1,:), p_ell(2,:), 'r')
	



