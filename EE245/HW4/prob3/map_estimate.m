function X_est  = map_estimate(X_est, z, dp_meas, descriptors, X_1_prior, P_prior, timesteps, Q, R)

J = [0 -1; 1 0];

% find the number of rows in A:
% the total number of landmark measurements
land_meas_count = 0;
for i = 1:length(z)
	land_meas_count = land_meas_count+length(z(i).descriptors);
end

% the total number of rows in A: 
rows_A = 3 + 3*(timesteps-1) + 2*land_meas_count

% maximum number of iterations
max_iter = 10;

for iter = 1:max_iter
	
	%initialize matrices
	A = zeros(rows_A,length(X_est));
	b = zeros(rows_A,1);
	
	% The prior 
	res = X_1_prior - X_est(1:3);
	% make sure the angle error is in the first cycle
	res(3) = atan2(sin(res(3)),cos(res(3)));

	% add the elements to A and b
	A(1:3,1:3) = sqrtm(inv(P_prior));
	b(1:3) = sqrtm(inv(P_prior))*res;
	
	% index to the next row of A where measurements should go
	row_ind = 4;
	
	% the odometry measurements: 
	% the square root inverse of Q
	sqrtQinv = inv(sqrtm(Q));
	for k = 1:timesteps-1
        
        % The current robot state and the next state
        X_k = X_est(3*k-2:3*k);
        X_k1 = X_est(3*k+1:3*k+3);
        phi_k = X_k(3);
        C = [cos(phi_k) -sin(phi_k);
		 sin(phi_k)  cos(phi_k)];
        
        % Find the residual
        res = zeros(3,1);
        res(1:2) = dp_meas(1:2,k) - C'*(X_k1(1:2) - X_k(1:2));
        res(3) = dp_meas(3,k) - (X_k1(3) - X_k(3));
        % make sure the angle error is in the first cycle
        res(3) = atan2(sin(res(3)),cos(res(3)));
        
        % Find the process Jacobian
        F_k = [ C' C'*J*(X_k1(1:2) - X_k(1:2)); 0 0 1];
        F_k1 = [ -C' zeros(2,1); 0 0 -1];
	  
		A(row_ind:row_ind+2,3*k-2:3*k+3) = [sqrtQinv*F_k sqrtQinv*F_k1];
		b(row_ind:row_ind+2) = - sqrtQinv*res;
		row_ind = row_ind+3;
	end
 
	% the square root inverse of R
	sqrtinvR = sqrtm(inv(R));
	% add the features measurements
	for i = 1:timesteps
        % The current robot state
        X_k = X_est(3*i-2:3*i);
        phi_k = X_k(3);
        C = [cos(phi_k) -sin(phi_k);
		 sin(phi_k)  cos(phi_k)];
     
		for j = 1:length(z(i).descriptors)
			
			% the index to the landmark being measured
			land_ind = find(descriptors == z(i).descriptors(j));
			
			%the corresponding columns in A
			A_land_ind = 3*timesteps + [2*land_ind-1:2*land_ind];
            
            % Find the residual
            res = z(i).measurements(:,j) - C'*(X_est(A_land_ind) - X_k(1:2));
			
			% find the measurement Jacobian
            H_r = [C' C'*J*(X_est(A_land_ind) - X_k(1:2))];
            H_l = -C';

			A(row_ind:row_ind+1, 3*i-2:3*i) = sqrtinvR*H_r;
			A(row_ind:row_ind+1, A_land_ind) = sqrtinvR*H_l;
			b(row_ind:row_ind+1) = - sqrtinvR*res;
			  
			row_ind = row_ind+2;
		end % end if meas_no>1
    end
 
	% solve the system
	A = sparse(A);
	sol = A\b;
	disp(sprintf('solution norm %e',norm(sol)))
    
%     figure, spy(A)
%     keyboard;
	
	
	% update the state estimate
	X_est = X_est + sol;
	
	if norm(sol)<1e-4
		break
	end
	
end
