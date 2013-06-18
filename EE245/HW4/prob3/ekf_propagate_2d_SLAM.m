function [xk1k ,Pk1k] = kf_propagate_2d(xkk, Pkk, dp_meas, dt, Q);

% find the number of landmarks that we have:
N = (length(xkk)-3)/2;

% propagate robot state estimates
C = [cos(xkk(3)) -sin(xkk(3));
	 sin(xkk(3))  cos(xkk(3))];
xk1k(1:3,1) = xkk(1:3) + [C*dp_meas(1:2);
	            		    dp_meas(3)];

% the landmarks don't move
xk1k(4:3+2*N) = xkk(4:3+2*N);

% compute the Jacobians for covariance propagation:
J = [0 -1; 
	 1 0];
Phi_R = [eye(2)  J*C*dp_meas(1:2);
	      0    0    1];
 
G_R = [C zeros(2,1);
	   0 0 1];

% covariance propagation
% the robot covariance:
Pk1k(1:3,1:3) = Phi_R*Pkk(1:3,1:3)*Phi_R' + G_R*Q*G_R';
% the cross terms:
Pk1k(1:3,4:3+2*N) = Phi_R*Pkk(1:3,4:3+2*N);
Pk1k(4:3+2*N,1:3) = Pkk(4:3+2*N,1:3)*Phi_R';
% the landmark covariance:
Pk1k(4:3+2*N,4:3+2*N) = Pkk(4:3+2*N,4:3+2*N);
