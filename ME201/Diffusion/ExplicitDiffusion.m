function Temp = ExplicitDiffusion(h, k, alpha, L, Ta, To, Tm_limit, dx, Fo);
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

if (Fo >= 0.5 )
    disp('Fo > 0.5. Change dt to get Fo = 0.1.');
    Fo = 0.1;    
end
dt=Fo*dx^2/alpha;

% number of points along the x axis
Nx = round(L/dx)+1;
% number of points along the t axis
Nt = round(Tm_limit/dt)+1;
% number of variables: Nx - two points at terminals
Nb = Nx - 2;

% prepare for propagation equation
e = Fo*ones(Nb-1,1);
f = (1-2*Fo)*ones(Nb,1);
g = Fo*ones(Nb-1,1);
% propagation matrix
Prop = diag(e,-1) + diag(f,0) + diag(g,1);

% adjust for boundary condition with Temp(0,all timestep)
Bi = h*dx/k;

% Iteratively solve for the wave equation
Temp = zeros(Nx,Nt);
% initialize for the first time step
Temp(:,1) = To-Ta;
for timestep = 1:Nt-1
    Temp(2:Nx-1, timestep+1) = Prop*Temp(2:Nx-1, timestep);
    Temp(2, timestep+1) = Temp(2, timestep+1) + Fo*Temp(1, timestep);
    Temp(Nx-1, timestep+1) = Temp(Nx-1, timestep+1) + Fo*(To-Ta);
    
    Temp(1,timestep+1) = Temp(2,timestep+1)/(1+Bi);
    Temp(Nx,timestep+1) = To-Ta;
end

Temp = Temp + Ta;

end