clear all;
close all;

% parameters
v0 = 1.0;

h=0.2;
t = 0:h:15;
% analytical solution
f = @(t,v) -3*t.*v./(1+t) + 2*(1+t).^3*exp(-t);
g = @(t) (1+t).^3.*exp(-t);

v1 = g(t);
% Explicit Euler method
[t_exeuler, v1_exeuler] = explicit_euler(f,[0 15],v0,length(t));
% Implicit Euler method
t_imeuler = t;
v1_imeuler = zeros(size(t_imeuler));
v1_imeuler(1) = v0;
for i=1:length(v1_imeuler)-1
    beta = 2*(1+t_imeuler(i+1)).^3*exp(-t_imeuler(i+1));
    alpha = 3*t_imeuler(i+1)./(1+t_imeuler(i+1));
    v1_imeuler(i+1) = (v1_imeuler(i)+beta*h)./(1+alpha*h);
end
% Trapezoidal method
t_trapz = t;
v1_trapz = zeros(size(t_trapz));
v1_trapz(1) = v0;
for i=1:length(v1_trapz)-1
    beta1 = 2*(1+t_trapz(i)).^3*exp(-t_trapz(i));
    alpha1 = 3*t_trapz(i+1)./(1+t_trapz(i+1));
    beta2 = 2*(1+t_trapz(i+1)).^3*exp(-t_trapz(i+1));
    alpha2 = 3*t_trapz(i+1)./(1+t_trapz(i+1));
    
    num = v1_trapz(i).*(1-alpha1*h/2) + (beta1+beta2)*h/2;
    dem = 1 + alpha2*h/2;
    v1_trapz(i+1) = num./dem;
end
% 2nd-order Runge-Kutta
[t_rk2, v1_rk2] = rk2(f,[0 15],v0,length(t));
% 4th-order Runge_Kutta
[t_rk4, v1_rk4] = rk4(f,[0 15],v0,length(t));

figure(11)
plot(t, v1,'b');
hold on
plot(t_exeuler, v1_exeuler,'g');
plot(t_imeuler, v1_imeuler,'g*');
plot(t_trapz, v1_trapz, 'c');
plot(t_rk2, v1_rk2,'r');
plot(t_rk4, v1_rk4,'r*');
legend( 'Analytical (Exact)', 'Explicit Euler', 'Implicit Euler', 'Trapezoidal', '2nd-order RK', '4th-order RK' );
xlabel('t');
ylabel('v');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=0.8;
t = 0:h:15;
v2 = g(t);
% Euler method
[t_exeuler, v2_exeuler] = explicit_euler(f,[0 15],v0,length(t));
% Implicit Euler method
t_imeuler = t;
v2_imeuler = zeros(size(t_imeuler));
v2_imeuler(1) = v0;
for i=1:length(v2_imeuler)-1
    beta = 2*(1+t_imeuler(i+1)).^3*exp(-t_imeuler(i+1));
    alpha = 3*t_imeuler(i+1)./(1+t_imeuler(i+1));
    v2_imeuler(i+1) = (v2_imeuler(i)+beta*h)./(1+alpha*h);
end
% Trapezoidal method
t_trapz = t;
v2_trapz = zeros(size(t_trapz));
v2_trapz(1) = v0;
for i=1:length(v2_trapz)-1
    beta1 = 2*(1+t_trapz(i)).^3*exp(-t_trapz(i));
    alpha1 = 3*t_trapz(i+1)./(1+t_trapz(i+1));
    beta2 = 2*(1+t_trapz(i+1)).^3*exp(-t_trapz(i+1));
    alpha2 = 3*t_trapz(i+1)./(1+t_trapz(i+1));
    
    num = v2_trapz(i).*(1-alpha1*h/2) + (beta1+beta2)*h/2;
    dem = 1 + alpha2*h/2;
    v2_trapz(i+1) = num./dem;
end
% 2nd-order Runge-Kutta
[t_rk2, v2_rk2] = rk2(f,[0 15],v0,length(t));
% 4th-order Runge_Kutta
[t_rk4, v2_rk4] = rk4(f,[0 15],v0,length(t));

figure(21)
plot(t, v2,'b');
hold on
plot(t_exeuler, v2_exeuler,'g');
plot(t_imeuler, v2_imeuler,'g*');
plot(t_trapz, v2_trapz, 'c');
plot(t_rk2, v2_rk2,'r');
plot(t_rk4, v2_rk4,'r*');
legend( 'Analytical (Exact)', 'Explicit Euler', 'Implicit Euler', 'Trapezoidal', '2nd-order RK', '4th-order RK' );
xlabel('t');
ylabel('v');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=1.1;
t = 0:h:15;
v3 = g(t);
% Euler method
[t_exeuler, v3_exeuler] = explicit_euler(f,[0 15],v0,length(t));
% Implicit Euler method
t_imeuler = t;
v3_imeuler = zeros(size(t_imeuler));
v3_imeuler(1) = v0;
for i=1:length(v3_imeuler)-1
    beta = 2*(1+t_imeuler(i+1)).^3*exp(-t_imeuler(i+1));
    alpha = 3*t_imeuler(i+1)./(1+t_imeuler(i+1));
    v3_imeuler(i+1) = (v3_imeuler(i)+beta*h)./(1+alpha*h);
end
% Trapezoidal method
t_trapz = t;
v3_trapz = zeros(size(t_trapz));
v3_trapz(1) = v0;
for i=1:length(v3_trapz)-1
    beta1 = 2*(1+t_trapz(i)).^3*exp(-t_trapz(i));
    alpha1 = 3*t_trapz(i+1)./(1+t_trapz(i+1));
    beta2 = 2*(1+t_trapz(i+1)).^3*exp(-t_trapz(i+1));
    alpha2 = 3*t_trapz(i+1)./(1+t_trapz(i+1));
    
    num = v3_trapz(i).*(1-alpha1*h/2) + (beta1+beta2)*h/2;
    dem = 1 + alpha2*h/2;
    v3_trapz(i+1) = num./dem;
end
% 2nd-order Runge-Kutta
[t_rk2, v3_rk2] = rk2(f,[0 15],v0,length(t));
% 4th-order Runge_Kutta
[t_rk4, v3_rk4] = rk4(f,[0 15],v0,length(t));

figure(31)
plot(t, v3,'b');
hold on
plot(t_exeuler, v3_exeuler,'g');
plot(t_imeuler, v3_imeuler,'g*');
plot(t_trapz, v3_trapz, 'c');
plot(t_rk2, v3_rk2,'r');
plot(t_rk4, v3_rk4,'r*');
legend( 'Analytical (Exact)', 'Explicit Euler', 'Implicit Euler', 'Trapezoidal', '2nd-order RK', '4th-order RK' );
axis([0 15 -1 10]);
xlabel('t');
ylabel('v');