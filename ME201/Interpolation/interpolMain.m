% x = 1900:10:1990;
% y = [75.995  91.972  105.711  123.203  131.669...
%      150.697  179.323  203.212  226.505  249.633];
% %[x' y']
% 
% xin = 1945;
% yout = interpolNeville(x,y,xin)
% yout = interpolCubicSpline(x,y,xin)

% %Compare with Matlab commands
% interp1(x,y, 1975,'spline')
% % or equivalently
% yout = spline(x,y,xin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Question 8a)
x = 1978:2:1992;
y = [12.0 12.7 13.0 15.2 18.2 19.8 24.1 28.1];

% Question 8b)
xb = [x(1:2) x(5:end)];
yb = [y(1:2) y(5:end)];

yout = interpolNeville(xb,yb,1980)
yout = interpolNeville(xb,yb,1982)

% Question 8c)
yout = interpolCubicSpline(xb,yb,1980)
yout = interpolCubicSpline(xb,yb,1982)

% Comment: not correct for 1982, but close.