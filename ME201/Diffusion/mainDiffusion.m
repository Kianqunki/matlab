% script to solve the wave equation
% Solving the wave equation dT/dt = alpha*d^2T/dx^2
% with boundary condition (BC): k*dT.dx = h(T-Ta); T = Ta at x = L;
% initial condition (IC): T(x,0) = To
% To, Ta: object and ambient temperature

clear all;
close all;

h = 0.1;
alpha = 0.1;
k = 0.1;
Ta = 300;
To = 3000;

L = 1; % meter
dx = 0.1;
Tm_limit = 10; % seconds

% dt = 1;
% % Explicit solution
% Temp1 = ExplicitDiffusion(h, k, alpha, L, Ta, To, Tm_limit, dx, dt);
% % Implicit solution
% Temp2 = ImplicitDiffusion(h, k, alpha, L, Ta, To, Tm_limit, dx, dt);
% % Analytical solution
% Temp3 = AnalyticalDiffusion(h, k, alpha, L, Ta, To, Tm_limit, dx, dt);

Fo = 0.1;
% Explicit solution
Temp1 = ExplicitDiffusion(h, k, alpha, L, Ta, To, Tm_limit, dx, Fo);
% Implicit solution
Temp2 = ImplicitDiffusion(h, k, alpha, L, Ta, To, Tm_limit, dx, Fo);
% Analytical solution
Temp3 = AnalyticalDiffusion(h, k, alpha, L, Ta, To, Tm_limit, dx, Fo);

% Plot
% define time moment to plot
Nprint = 2;
x = (0:size(Temp1,1)-1)*dx;

h=figure(3);
plot(x,Temp1(:,Nprint),'-r');
hold on;
plot(x,Temp2(:,Nprint),'-g');
plot(x,Temp3(:,Nprint),'-b');
grid on;
xlabel('Distance'); ylabel('Temperature');
legend('Explicit', 'Implicit', 'Analytical',2);

% Diffusion(h,k,alpha,L,Tm_limit,(Nprint-1)/Tm_limit,Ta,To);
