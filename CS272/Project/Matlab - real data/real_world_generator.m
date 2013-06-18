function [X, v_m, omega_m, dt, z, X_L] = real_world_generator()

load('data_set.mat');
load('beac_juan3.mat');

omega_m = Steering;

% num=length(Time);
% nextTime=zeros(1,num);
% nextTime(1:num-1) = Time(2:num);
% nextTime(num) = Time(num);
% nextDiff = abs(nextTime-Time);
% idx=find(next>0);
% min(next(idx)) % => give 0.001 seconds
% idx=find(next<1);
% max(next(idx))  %=> give 0.103 seconds

num=length(Time_VS);
nextTime=zeros(1,num);
nextTime(1:num-1) = Time_VS(2:num);
nextTime(num) = Time_VS(num);
dt = nextTime-Time_VS;
v_m = Velocity;

% find true position from GPS
% plotting will be different
X = [GPSLat; GPSLon];   %Latitude = x, Longtitude = y

% find the measurements from laser
lastIdx = 1;
for k = 1:length(TimeLaser)
    t = TimeLaser(k);
    diff_t = abs(t - Time_VS(:));
    [t_min, idx] = min(diff_t);
    
    for i=lastIdx:idx-1
        z(i).measurements=[];
    end
    
    % read in the measurements from laser sensor
    laserIdx = find(Intensity(k,:)>0);
       
    z(idx).measurements=zeros(2,length(laserIdx));
    alpha=laserIdx/2;   %index from 0..361 for 0..180
    for i=1:length(laserIdx)
        distance=Laser(k,laserIdx(i));
        z(idx).measurements(1,i)=distance*cosd(alpha(i));
        z(idx).measurements(2,i)=distance*sind(alpha(i));
    end

    lastIdx=idx+1;
end

X_L = estbeac';  % true position of the beacons.
end