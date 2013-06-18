function Tanal = AnalyticalDiffusion(h, k, alpha, L, Ta, To, Tm_limit, dx, Fo);
%Solving the wave equation dT/dt = alpha*d^2T/dx^2
% with boundary condition (BC): k*dT.dx = h(T-Ta); T = Ta at x = L;
% initial condition (IC): T(x,0) = To
% To, Ta: object and ambient temperature

% Input: 
% h,k,alpha,L,Ta,To: equation parameters
% Tm_limit: length of time for simulation
% dx, dt: discretization parameters
% Output:
% Temp: temperature values in 2D array

dt=Fo*dx^2/alpha;

% number of points along the x axis
Nx = round(L/dx)+1;
% number of points along the t axis
Nt = round(Tm_limit/dt)+1;
% number of variables: Nx - two points at terminals
Nb = Nx - 2;

x = (0:Nx-1)*dx;
tm = (0:Nt-1)*dt;

Tanal = zeros(Nx,Nt);
for it=1:Nt
    for ix=1:Nx
        aa=h/k; bb=sqrt(alpha*tm(it));
        cc=x(ix)/(2*bb+1.0e-08);
        Tanal(ix,it)=erf(cc)+...
            exp(aa*x(ix)+(aa*bb)^2)*(1-erf(cc+aa*bb));
    end
end
Tanal=(To-Ta)*Tanal+Ta;