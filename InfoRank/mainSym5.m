% Case 5: 4 robot poses, 2 landmarks, always observed, with odometry
close all
clear all

% define components of robot poses at time k: \hat{x}_r(k)
% suffix m (marginalized) for time k, suffix a (appended) for time k'
syms x1m y1m phi1m real;
syms x2m y2m phi2m real;
syms x3m y3m phi3m real;
syms x4m y4m phi4m real;
% define components of robot poses at time k': \hat{x}_r(k')
syms x1a y1a phi1a real;
syms x2a y2a phi2a real;
syms x3a y3a phi3a real;
syms x4a y4a phi4a real;

% define components of landmarks at time k: \hat{x}_l(k) = \hat{p}_l(k)
syms m1m n1m real;
syms m2m n2m real;
% define components of landmarks at time k': \hat{x}_l(k') = \hat{p}_l(k')
syms m1a n1a real;
syms m2a n2a real;

% define the robot pose vectors
r_m = [x1m y1m phi1m; x2m y2m phi2m; x3m y3m phi3m; 0 0 0]';
r_a = [ 0 0 0; x2a y2a phi2a; x3a y3a phi3a; x4a y4a phi4a]';
r_all = [ x1a y1a phi1a; x2a y2a phi2a; x3a y3a phi3a; x4a y4a phi4a]';
% define the landmark vectors
l_m = [m1m n1m; m2m n2m]';
l_a = [m1a n1a; m2a n2a]';
l_all = [m1a n1a; m2a n2a]';

% define S_m(k)
S_m = [1 1; 1 1; 1 1; 0 0];
% define S_a(k')
S_a = [0 0; 1 1; 1 1; 1 1];
% define union of S_m(k) and S_a(k')
S = double(S_m | S_a);

% number of components for each landmarks
params.lmDim = 2;
params.rpDim = 3;
params.J = [0 -1; 1 0];
% number of robot poses
params.K = size(S,1);
% number of landmarks
params.N = size(S,2);

% Compute the rank of A**
Add = sym(zeros(params.rpDim*params.K+params.lmDim*params.N));
% Add the measurement Jacobian
for i=1:size(S,1)
    for j=1:size(S,2)    
        H = symMeasJacobian(i,j,S,r_all,l_all,params);
        Add = Add + simplify(H'*H);
    end
end
% Add the odemetry propagation Jacobian
for i=1:size(S,1)
    H = symPropJacobian(i, S, r_all, params);
    H
    keyboard;
    Add = Add + simplify(H'*H);
end
Add = simplify(Add);

% % Compute the rank of A*
% Ad = sym(zeros(params.rpDim*params.K+params.lmDim*params.N));
% for i=1:size(S_m,1)
%     for j=1:size(S_m,2)    
%         H = symMeasJacobian(i,j,S_m,r_m,l_m,params);
%         Ad = Ad + simplify(H'*H);
%     end
% end
% for i=1:size(S_a,1)
%     for j=1:size(S_a,2)    
%         H = symMeasJacobian(i,j,S_a,r_a,l_a,params);
%         Ad = Ad + simplify(H'*H);
%     end
% end
% Ad = simplify(Ad);
% null(Ad)
% rref(Ad)