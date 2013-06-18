clear all;
close all;

% data
x=[1.1,1.65,3.05,4.8,6.9,9.6,13.1]';
t=[0,0.5,1.0,1.5,2.0,2.5,3.0]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% batch LS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
std = 0.1;
% build matrix W
W = (1/std^2)*eye(size(x,1));
% build matrix H
t0 = ones(size(x));
t1 = t;
t2 = t.^2;
H = [t0 t1 t2];

disp( 'Batch LS Estimation:');
theta1 = inv(H'*W*H)*H'*W*x

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RLS in Information form
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p=3; % number of parameters

% initialize by batch LS with partial data
W_n = (1/std^2)*eye(p);
H_n = H(1:p,:);
theta_n = inv(H_n'*W_n*H_n)*H_n'*W_n*x(1:p);
I_n = H_n'*W_n*H_n;

for i=p+1:size(x,1)
    % just for clarity purpose
    h_n = H(i,:)';
    theta_prev = theta_n;
    I_prev = I_n;
    
    % information matrix update
    I_n = I_prev + h_n*h_n'/std^2;
    % invert the matrix
    P_n = inv(I_n);
    % compute K_n
    K_n = P_n*h_n/std^2;
    % estimator update
    theta_n = theta_prev + K_n*(x(i) - h_n'*theta_prev);
end

disp( 'RLS in Information form:');
theta_n 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RLS in Covariance form
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear theta_n % NO CHEATING
p=3; % number of parameters

% initialize by batch LS with partial data
W_n = (1/std^2)*eye(p);
H_n = H(1:p,:);
theta_n = inv(H_n'*W_n*H_n)*H_n'*W_n*x(1:p);
Sigma_n = inv(H_n'*W_n*H_n);

for i=p+1:size(x,1)
    % just for clarity purpose
    h_n = H(i,:)';
    theta_prev = theta_n;
    Sigma_prev = Sigma_n;
    
    % compute K_n, equation 8.47
    K_n = Sigma_prev*h_n/(std^2 + h_n'*Sigma_prev*h_n);
    % estimator update, equation 8.46
    theta_n = theta_prev + K_n*(x(i) - h_n'*theta_prev);
    % covariance update, equation 8.48
    Sigma_n = (eye(p) - K_n*h_n')*Sigma_prev;
end

disp( 'RLS in Covariance form:');
theta_n

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Batch LS Estimation:
% 
% theta1 =
% 
%     1.0762
%     0.8464
%     1.0452
% 
% RLS in Information form:
% 
% theta_n =
% 
%     1.0762
%     0.8464
%     1.0452
% 
% RLS in Covariance form:
% 
% theta_n =
% 
%     1.0762
%     0.8464
%     1.0452