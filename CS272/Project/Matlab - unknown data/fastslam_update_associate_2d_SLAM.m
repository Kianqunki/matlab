function [particles, descriptors] = fastslam_update_associate_2d_SLAM(particles, z, R, descriptors,X_L)
% fastslam_update_2d_SLAM: the two update steps in FastSLAM v1.0
% namely, Update landmark estimates, Calculate importance weight

% z.descriptors(i) is still generated but not used in this function. Only
% used for debugging and evaluation.

particleNum = length(particles);
J = [0 -1; 1 0];
T1=chi2inv(0.95,2);
T2=chi2inv(0.99999,2);

for i=1:size(z.measurements,2)
    for m=1:particleNum         %number of particles
        landmarkNum = size(particles(m).Xl,2); %number of landmarks
        check = zeros(1, landmarkNum );
        
        pose = particles(m).Xr;
        C = [cos(pose(3)) -sin(pose(3));
            sin(pose(3)) cos(pose(3))];
        
        for n=1:length(check)   %number of landmarks
            H = C';         % corresponding to G in the thesis
            z_hat = C'*(particles(m).Xl(:,n) - pose(1:2));
            r = z.measurements(:,i) - z_hat;   % measurement innovation
            S = H*particles(m).Pl(:,:,n)*H'+R;
            check(n) = r'*inv(S)*r;
        end
        
        idx=find(check>T2);
        if (length(idx)==length(check))
            % out of bounds for all known landmarks => new landmark            
            
            % initialize the fields
            if (landmarkNum == 0)
                particles(m).Xl = zeros(2,1);
                particles(m).Pl = zeros(2,2,1);
            end
            
            % add estimated position of the new landmark
            particles(m).Xl(:,landmarkNum+1) = pose(1:2) + C*z.measurements(:,i);

            % add covariance for the new landmark.
            % section 3.4.3 in the thesis
            % find Jacobian for measurment model
            H = C';
            % if H is rotation matrix, then inv(H'*inv(R)*H) = H'*R*H;
            particles(m).Pl(:,:,landmarkNum+1) = H'*R*H;
            
            % weight update is not needed because weights are unchanged
            % DEBUG
            particles(m).dict = [particles(m).dict z.descriptors(i)];
        end
        
        idx=find(check<T1);
        if (length(idx)==1)
            % in the boundary of a single known landmark
            z_hat = C'*(particles(m).Xl(:,idx) - pose(1:2));
            r = z.measurements(:,i) - z_hat;   % measurement innovation
            H = C';
            S = H*particles(m).Pl(:,:,idx)*H'+R;
            K = particles(m).Pl(:,:,idx)*H'/S;
            particles(m).Xl(:,idx) = particles(m).Xl(:,idx) + K*r;
            particles(m).Pl(:,:,idx) = (eye(size(H)) - K*H)*particles(m).Pl(:,:,idx);
            
%             % try to moderate the landmark error
%             T = 0.1;
%             if (min(diag(particles(m).Pl(:,:,idx))) < T )
%                 particles(m).Pl(:,:,idx) = particles(m).Pl(:,:,idx).*T/min(diag(particles(m).Pl(:,:,idx)));
%             end
            
            % weight update
            % Calculate the importance weight
            % section 3.3.3 in the thesis, equation 3.37
            % NOTE: the equation 3.37 is slightly wrong or different way of
            % notations
            gauss = exp(-0.5*r'*inv(S)*r)/2*pi*sqrt(det(S));
            particles(m).weight = particles(m).weight*gauss;
        end
        
        % otherwise, do nothing: no association or no new landmark        
    end
end


