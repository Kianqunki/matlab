function [t,y] = rk4(f,tspan,y0,N)
% Fourth-order Runge-Kutta

m = length(y0);
t = linspace(tspan(1),tspan(2),N+1);
y = zeros(m,N+1);
h = (tspan(2)-tspan(1))/N;
y(:,1) = y0;

for i = 1:N
    k1 = h*f(t(i),y(:,i));
    k2 = h*f(t(i)+h/2,y(:,i)+k1/2);
    k3 = h*f(t(i)+h/2,y(:,i)+k2/2);
    k4 = h*f(t(i)+h/2,y(:,i)+k3);
    y(:,i+1) = y(:,i) + k1/6 + (k2+k3)/3 + k4/6;
end
