function [particles, descriptors] = fastslam_update_2d_SLAM(particles, z, R, descriptors,X_L)
% fastslam_update_2d_SLAM: the two update steps in FastSLAM v1.0
% namely, Update landmark estimates, Calculate importance weight

particleNum = length(particles);
%landmarkNum = length(descriptors);
landmarkNum = size(particles(1).Xl,2);
J = [0 -1; 1 0];

for j=1:length(z.descriptors)
    idx = find(descriptors==z.descriptors(j),1);
    %the second parameter 1 is for faster processing

    % if the landmark is new
    if (isempty(idx))
        for i=1:particleNum
            % initialize the fields
            if (landmarkNum == 0)
                particles(i).Xl = zeros(2,1);
                particles(i).Pl = zeros(2,2,1);
            end

            % add the new features
            pose = particles(i).Xr;
            C = [cos(pose(3)) -sin(pose(3));
                sin(pose(3)) cos(pose(3))];

            % add estimated position of the new landmark
            particles(i).Xl(:,landmarkNum+1) = pose(1:2) + C*z.measurements(:,j);

            % add covariance for the new landmark.
            % section 3.4.3 in the thesis
            % find Jacobian for measurment model
            H = C';
            % if H is rotation matrix, then inv(H'*inv(R)*H) = H'*R*H;
            particles(i).Pl(:,:,landmarkNum+1) = H'*R*H;
            %particles(i).Pl(:,:,landmarkNum+1) = 0.02*eye(2);
        end
        
        % update the list of landmark IDs
        descriptors = [descriptors z.descriptors(j)];
        landmarkNum = landmarkNum+1;
        
%         %DEBUG
%         X_L(:,z.descriptors(j))
%         particles(1).Pl(:,:,landmarkNum)
%         input('Wait here');
    else % if the landmark is already discovered
        for i=1:particleNum
            % Udpate the EKFs
            % section 3.3.2 in the thesis
            pose = particles(i).Xr;
            C = [cos(pose(3)) -sin(pose(3));
                sin(pose(3)) cos(pose(3))];
            
            z_hat = C'*(particles(i).Xl(:,idx) - pose(1:2));
            r = z.measurements(:,j) - z_hat;   % measurement innovation
            H = C';
            S = H*particles(i).Pl(:,:,idx)*H'+R;
            K = particles(i).Pl(:,:,idx)*H'/S;
            particles(i).Xl(:,idx) = particles(i).Xl(:,idx) + K*r;
            particles(i).Pl(:,:,idx) = (eye(size(H)) - K*H)*particles(i).Pl(:,:,idx);
            
            % Calculate the importance weight
            % section 3.3.3 in the thesis, equation 3.37
            % NOTE: the equation 3.37 is slightly wrong.
            gauss = exp(-0.5*r'*inv(S)*r)/2*pi*sqrt(det(S));
            particles(i).weight = particles(i).weight*gauss;
        end
        
%         %DEBUG
%         particles(i).Pl(:,:,idx)
%         input('Wait here');
    end
end


