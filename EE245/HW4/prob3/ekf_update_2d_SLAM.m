function [X_est,P_est, descriptors]= kf_update_2d(X_est, P_est, z, R, descriptors)

% this update function uses sequential updates 

% see how many measurements we have:
M = length(z.descriptors);

for i = 1:M
	
	% check if this landmark is already in the state vector:
	ind = find (z.descriptors(i) == descriptors);
	
	if length(ind) == 1
		% then this landmark is already in the state, carry out an update 
		
		% the estimated rotation matrix:
		C_est = [cos(X_est(3))    -sin(X_est(3));
			     sin(X_est(3))     cos(X_est(3))];
		
		J = [0 -1; 1 0];
		 
		% the estimate for this landmark:
		X_L_est = X_est(3+2*ind-1:3+2*ind);
		
		% build the measurement Jacobian:
		H_R = C_est'* [-eye(2)  -J*(X_L_est - X_est(1:2))];
		H_L = C_est';

		H = zeros(2,length(X_est));
		H(:,1:3) = H_R;
		H(:,3+2*ind-1:3+2*ind)= H_L;

		% carry out the update:
		
		% the expected measurement
		z_hat = C_est'*(X_L_est-X_est(1:2));
				
		% the residual covariance matrix:
		S = H*P_est*H' + R;

		% the Kalman gain
		K = P_est*H'*inv(S);
 
		% the measurement residual
		r = z.measurements(:,i) - z_hat;

		% update the state estimate
		X_est = X_est + K*r;

		% update the covariance matrix
		P_est = (eye(length(P_est)) - K*H)*P_est*(eye(length(P_est)) - K*H)' + K*R*K';

	else
		
		%initialize a new landmark (following class notes)
		descriptors = [descriptors z.descriptors(i)];
 		
		% the estimated rotation matrix:
		C_est = [cos(X_est(3))    -sin(X_est(3));
			     sin(X_est(3))     cos(X_est(3))];
 
		% the landmarks estimate
   	    X_L_est = X_est(1:2) + C_est*z.measurements(:,i);

		% add it in the state vector:
		X_est = [X_est; X_L_est];
		
		% now we need to augment the covariance matrix
	
		% build the measurement Jacobian:
		J = [0 -1; 1 0];
		H_R = C_est'* [-eye(2)  -J*(X_L_est - X_est(1:2))];
		H_L = C_est';

		% the 2x2 covariance matrix of the new landmark:
		P_est(end+1:end+2,end+1:end+2) = H_L'*(H_R*P_est(1:3,1:3)*H_R'+R)*H_L;
		
		% the cross terms:
		P_est(1:end-2,end-1:end) = -P_est(1:end-2,1:3) * H_R'*H_L;
		P_est(end-1:end,1:end-2) = P_est(1:end-2,end-1:end)';
		
		
	end
	
	
end

 
