function dp = scan_match( scan1, scan2 )
% ICP algorithm for motion estimation

% figure(1)
% plot(scan1.points(2,:),scan1.points(1,:),'*');
% hold on
% axis([-3 3 0 10])
% figure(2)
% plot(scan2.points(2,:),scan2.points(1,:),'*');
% axis([-3 3 0 10])

% The relative motion between the two scans
dp = scan2.dp_est;
npts2 = size(scan2.points,2);
npts1 = size(scan1.points,2);
projected_scan1 = zeros(2,npts2);

for iter=1:6
    % Find the projection of scan 2 onto scan 1
    phi = dp(3);
    C = [cos(phi)  -sin(phi);
        sin(phi)    cos(phi)];
    for j=1:npts2
        projected_scan1(:,j) = dp(1:2) + C*scan2.points(:,j);
    end
%     figure(1)
%     plot(projected_scan1(2,:),projected_scan1(1,:),'r*');
%     hold off
    
    % For each point in the scan 2 (projected on scan 1)
    min_dist_2 = zeros(1,npts2);
    % Index for points IN scan 1 for each point in scan2
    match_index_2 = zeros(1,npts2);
    for j=1:npts2
        dist = zeros(1,npts1);
        for k=1:npts1
            diff = projected_scan1(:,j) - scan1.points(:,k);
            dist(k) = norm(diff);
        end
        
        % find the closest point to that point
        [min_dist_2(j) match_index_2(j)] = min(dist);
    end
    
    % for each point in the scan 1
    min_dist_1 = zeros(1,npts1);
    % Index for points in scan 2 for each point in scan1
    match_index_1 = zeros(1,npts1);
    for j=1:npts1
        dist = zeros(1,npts2);
        for k=1:npts2
            diff = projected_scan1(:,k) - scan1.points(:,j);
            dist(k) = norm(diff);
        end
        
        % find the closest point to that point
        [min_dist_1(j) match_index_1(j)] = min(dist);
    end
    
    % find the mutual matches, and save the matching index into match_index_2
    for j=1:npts2
        % if there is mutual matches and the distance is less than 0.5 meter
        if ( match_index_1( match_index_2(j) ) == j && min_dist_2(j) < 0.5 )
            % we find a match, keep it in match_index_2
            
            %         % VERIFY: distance must be the same
            %         min_dist_1(match_index_2(j))
            %         min_dist_2(j)
            %         figure(10)
            %         plot(projected_scan1(2,j),projected_scan1(1,j),'r*');
            %         hold on
            %         plot(scan1.points(2,match_index_2(j)),scan1.points(1,match_index_2(j)),'*');
            %         axis([-3 3 0 10])
            %         pause
        else
            % not matched, put negative flag
            match_index_2(j) = -1;
        end
    end
    
    % The correspondences will be 'index' and match_index_2(index)
    index = find( match_index_2 > 0 );
    
    last_dp = dp;
	% Find ML estimate based on the found correspondences
    dp = estimate_pose_GN( scan1.points(:,match_index_2(index)), scan2.points(:,index), scan1.Cov(:,:,match_index_2(index)), scan2.Cov(:,:,index), last_dp, length(index) );
    
    if ( norm(dp-last_dp) < 1e-3 )
        break
    end
end

% close all

end