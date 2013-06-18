function [t,y] = rk2(f,tspan,y0,N)
% Second-order Runge-Kutta

m = length(y0);
t = linspace(tspan(1),tspan(2),N+1);
y = zeros(m,N+1);
h = (tspan(2)-tspan(1))/N;
y(:,1) = y0;

alpha = 1/2;

for i = 1:N
    k1 = h*f(t(i),y(:,i));
    k2 = h*f(t(i)+alpha*h,y(:,i)+alpha*k1);
    y(:,i+1) = y(:,i) + (1-1/2/alpha)*k1 + k2/2/alpha;
end
