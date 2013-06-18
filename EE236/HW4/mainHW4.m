clear all
close all
randn('seed',1);

% KNOWN INFO
% define the antennas
A(:,1) = [0, 0]';
A(:,2) = [100, 0]';
A(:,3) = [0, 100]';
% signal propagation speed
c = 1000;

% UNKNOWN INFO
% define the source
Ps = [25, 60]';
% define the noise
noise_std = 0.0001;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate data: arrival time t
num = size(A,2); % number of data
t = zeros(1,num);
noise = noise_std*randn(1,num);
for i=1:num
    rel_dist = A(:,i) - Ps;
    t(:,i) = norm(rel_dist)/c + noise(:,i);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% process measurements
% initialize the source position
P = [1,0]';
deltaP = [1000, 1000]';

convCond = 1e-10; % convergence condition
maxIteration = 50; % maximum number of iterations

for i=1:maxIteration
    [deltaP, covP] = estimate_BLUE( A, c, noise_std, t, P );
    P = P + deltaP;
    
    deltaP
    if ( norm(deltaP) < convCond )
        % print the final covariance matrix at convergence point
        covP
        P
        break;
    end
end