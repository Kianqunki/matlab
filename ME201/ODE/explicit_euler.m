function [t,y] = explicit_euler(f,tspan,y0,N)
% Explicit Euler

m = length(y0);
t = linspace(tspan(1),tspan(2),N+1);
y = zeros(m,N+1);
h = (tspan(2)-tspan(1))/N;
y(:,1) = y0;

for i = 1:N
    y(:,i+1) = y(:,i) + h*f(t(i),y(:,i));
end