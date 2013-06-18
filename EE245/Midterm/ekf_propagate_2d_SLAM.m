function [xk1k ,Pk1k] = kf_propagate_2d(xkk, Pkk, vm, wmega, dt, sigma_v, sigma_omega);

% find the number of landmarks that we have:
N = (length(xkk)-3)/2;

% propagate robot state estimates
xk1k(1:3,1) = xkk(1:3) + [vm*dt*cos(xkk(3));
	             vm*dt*sin(xkk(3));
				 wmega*dt];
			 
% the landmarks don't move
xk1k(4:3+2*N) = xkk(4:3+2*N);

% compute the Jacobians for covariance propagation:
Phi_R = [1  0  -vm*dt*sin(xkk(3));  
	     0  1   vm*dt*cos(xkk(3));  
	     0  0     1];

G_R = [-dt*cos(xkk(3)) 0;
     -dt*sin(xkk(3)) 0;
     0             -dt];

% the covariance matrix of the process noise
Q = [sigma_v^2     0 ;
      0    sigma_omega^2];

% covariance propagation
% the robot covariance:
Pk1k(1:3,1:3) = Phi_R*Pkk(1:3,1:3)*Phi_R' + G_R*Q*G_R';
% the cross terms:
Pk1k(1:3,4:3+2*N) = Phi_R*Pkk(1:3,4:3+2*N);
Pk1k(4:3+2*N,1:3) = Pkk(4:3+2*N,1:3)*Phi_R';
% the landmark covariance:
Pk1k(4:3+2*N,4:3+2*N) = Pkk(4:3+2*N,4:3+2*N);
