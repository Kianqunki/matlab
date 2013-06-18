clear all
close all

load scans
 
for i = 1:length(scan_struct)-1
	i % printout to make sure things are moving
	dp_est = scan_match(scan_struct(i), scan_struct(i+1));
	% save the estimate from ICP in the structure for plotting later
	scan_struct(i+1).dp_est_icp = dp_est;
end


%%%% PLOT THE RESULTS

% put all scans in one map
% initialize
all_points = [];
robot_trajectory = zeros(2,1);
% the robot pose in the global frame
g_est = [0 0 0]';
for i = 2:length(scan_struct) 
	% find the current global pose
	phi = g_est(3);
	dp = scan_struct(i).dp_est_icp(1:2);
	dphi = scan_struct(i).dp_est_icp(3);
	C = [cos(phi)  -sin(phi);
	 	 sin(phi)    cos(phi)];
	g_est = g_est + [C*dp;dphi];
	% save position for plotting
	robot_trajectory = [robot_trajectory g_est(1:2)];
	
	% transform the scan points to the global frame using this estimate
	phi = g_est(3);
	C = [cos(phi)  -sin(phi);
	 	 sin(phi)    cos(phi)];
	p  = scan_struct(i).points;
	p_in_g = diag(g_est(1:2))*ones(size(p))+C*p;
	% save for plotting
	all_points = [all_points p_in_g];
end
  

% put all scans in one map, using the odometry estimates
% initialize
all_points_od = [];
robot_trajectory_od = zeros(2,1);
	
g_est = [0 0 0]';
for i = 2:length(scan_struct) 
	% find the current global  pose
	phi = g_est(3);
	dp = scan_struct(i).dp_est(1:2);
	dphi = scan_struct(i).dp_est(3);
	C = [cos(phi)  -sin(phi);
	 	 sin(phi)    cos(phi)];
	
	g_est = g_est + [C*dp;dphi];
	robot_trajectory_od = [robot_trajectory_od g_est(1:2)];

	% transform the scan points to the global frame using this estimate
	phi = g_est(3);
	C = [cos(phi)  -sin(phi);
	 	 sin(phi)    cos(phi)];
	p  = scan_struct(i).points;
	p_in_g = diag(g_est(1:2))*ones(size(p))+C*p;
	% save for plotting
	all_points_od = [all_points_od p_in_g];

end
  
%plot the points and the trajectories
figure
plot(all_points(1,:),all_points(2,:),'*')
hold on
plot(all_points_od(1,:),all_points_od(2,:),'r*')
plot(robot_trajectory(1,:),robot_trajectory(2,:))
plot(robot_trajectory_od(1,:),robot_trajectory_od(2,:),'r')
axis equal
legend('ICP', 'Odometry only', 'ICP trajectory', 'Odometry trajectory')