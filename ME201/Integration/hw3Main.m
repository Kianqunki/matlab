clear all;
close all;
format long;

% starting point of the line source
X = 0;
Y1 = 0;
% length of the line source
L = 1000;
% wind angle in degree
theta = 0:1:90;
% emission rate of the line source
q = 50; % don't know U, use this to represent q/U.
% position of the receptor
Xr = 1000;
Yr = 500;
% for minimum of 10^5 point sources
order = ceil(log2(10^5));

% plot the outputs
outputRomberg = zeros(length(theta),1);
outputGaussian = zeros(length(theta),1);
outputHV = zeros(length(theta),1);

for i=1:length(theta)
    outputRomberg(i) = hw3Romberg(X,Y1,L,Xr,Yr,theta(i),q,order);
    outputGaussian(i) = hw3Gaussian(X,Y1,L,Xr,Yr,theta(i),q,order);
end
outputHV = hw3AnalyticalHV(X,Y1,L,Xr,Yr,theta,q);
outputLE = hw3AnalyticalLE(X,Y1,L,Xr,Yr,theta,q);

figure(2)
norm = outputRomberg(1);
semilogy(theta, outputRomberg./norm, 'bx');
hold on;
norm = outputGaussian(1);
plot(theta, outputGaussian./norm, 'g+');
norm = outputLE(1);
semilogy(theta, outputLE./norm, 'r*');
norm = outputHV(1);
semilogy(theta, outputHV./norm, 'cs');
legend( 'Numerical Romberge', 'Numerical Gaussian', 'Analytical LE', 'Analytical HV', 'Location', 'SouthWest' );
%legend( 'Numerical Romberge', 'Analytical LE', 'Analytical HV', 'Location', 'SouthWest' );

% add axis labels
axis([0  90  1e-4  1.1])
xlabel ('Wind Angle (degrees)')
ylabel ('Normalized Concentration')

% % compare arbitrary numerical output
% theta = 30;
% % Romberg integration
% outputRomberg = hw3Romberg(X,Y1,L,Xr,Yr,theta,q,order)
% % Gaussian quadrature
% outGaussian = hw3Gaussian(X,Y1,L,Xr,Yr,theta,q,order)
% %compare with Matlab function
% outMatlab = quadl(@(Y)hw3SourceContrib(X,Y,Xr,Yr,theta,q),Y1,Y1+L)
