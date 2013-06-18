clear all
close all

seedNum = 1;
randn('seed',seedNum)
rand('seed',seedNum)

% define parameters
% true parameters
params.trueTheta = [1.0, -0.5, 0.25, -0.125, 0.0625]';
% standard deviation of Gaussian noise
params.sigma = 0.1;
% dimension of parameters
params.p = 5;

% number of observations
N = 150;

% generate the observations and input
[U, X] = data_generator( params, N );

% initialize before processing
% structure to save the processing results
thetaErr_sav = zeros(params.p,length(X));
C_sav = zeros( params.p, params.p, length(X) );

% processing
for m=5:length(X)
    [theta, C] = estimate_linear(X(1:m), U(1:m), params);
    
    % save for plotting
    thetaErr_sav(:,m) = params.trueTheta - theta;
    C_sav(:,:,m) = C;
end

% plotting
std = zeros(params.p,N);
for i=1:params.p
    % find the CRLB
    std_theta(i,:) = sqrt(C_sav(i,i,:));
    
    % plot the error with +-sqrt of CRLB
    figure(i)
    p1 = plot( thetaErr_sav(i,:) );
    hold on
    p2 = plot( 3*std_theta(i,:), 'red' );
    plot(-3*std_theta(i,:), 'red' );
    
    legend([p1(1) p2(1)], 'theta_i error', '\pm 3 sqrt CRLB')
end