function [C, x, z] = hw5PDEsolver(u, Kz , L, u_star, params)
% output x, z for index searching

% get the parameters
z0 = params.z0;
zm = params.zm;
zmax = params.zmax;
h = params.h;
xmax = params.xmax;

% number of steps in x direction
Nx = 200;
dx = xmax/Nx;
x = 0:dx:xmax;

z = [z0 logspace(log10(1.01*z0),log10(zmax),500)];
Nz = length(z);

uVal = zeros(Nz,1);
KzVal = zeros(Nz,1);

% compute all u and Kz
for j=1:Nz
    uVal(j) = u(z(j), z0, L, u_star);
    KzVal(j) = Kz(z(j), z0, L, u_star);
end

% intial value conditions
C = zeros(Nx,Nz);
% uh = u(h, z0, L, u_star);
% % approximate the initial condition with a Gaussian
% coeff = sum(fspecial('gaussian',3));
% [diff zIdx] = min(abs(h-z));
% C(1,zIdx-1:zIdx+1) = coeff(:)./uh;

% [diff hIdx] = min(abs(h-z));
% Q = (z(hIdx+1)-z(hIdx))*(uVal(hIdx+1)+uVal(hIdx))/2;
Q = 1;

ubar = (uVal(1)+uVal(2))/2;
C(1,1) = Q/ubar/(z(2)-z(1));
C(1,2) = C(1,1);

for i=2:Nx
    % find the vertical distribution
    f = zeros(Nz,1);
    g = zeros(Nz,1);
    e = zeros(Nz,1);
    b = zeros(Nz,1);
    
    for j=2:Nz-1
        % KzMinus = (KzVal(j-1)+KzVal(j))/2;
        num = (KzVal(j-1)+KzVal(j))*dx/2;
        dem = (z(j+1)-z(j))*(z(j+1)-z(j-1))*uVal(j);
        alphaMinus = num/dem;
        % KzPlus = (KzVal(j)+KzVal(j+1))/2;
        num = (KzVal(j)+KzVal(j+1))*dx/2;
        dem = (z(j)-z(j-1))*(z(j+1)-z(j-1))*uVal(j);
        alphaPlus = num/dem;
        
        e(j) = -alphaMinus;
        f(j) = 1+alphaPlus+alphaMinus;
        g(j) = -alphaPlus;
        b(j) = C(i-1,j);
    end
    
    % boundary conditions
    f(1) = 1;
    g(1) = -1;
    f(Nz) = 1;
%     e(Nz) = -1;
    
    % solve for C(i,:)
    C(i,:) = Tridiag_Solver(e,f,g,b);
%     C(i,Nz) = C(i,Nz-1);
%     C(i,1) = C(i,2); %/(1+beta)
end

end