% A Numerical Study on the dispersion of Passive contaminents from a
% continous source in the atmospheric surface layer by A.P.Van Ulden (1978)
function output=surface_release(Q,hs,u,k,Lx,Lz,z0,L,ustar)

% Q is the source strength in  mg/s
% u is the function  of velocity 
% k is the eddy diffusivity in m^2/s 
% Lx and Lz are the limits of the concentration field

delx=1; % X steps 

x=[0:delx:Lx]; % Downwind Distance Grid points

Nx=length(x);

 z1=hs-logspace(-1,log10(hs-z0),1e3);
 z2=logspace(log10(hs),log10(Lz),1e3);
 z=[z1(length(z1):-1:1) z2];

Nz=length(z);

delz=z(2:Nz)-z(1:Nz-1);

N=Nz-1;

for j=1:N-2
    
    K(j,1)=(k(ustar,L,z(j+1))+k(ustar,L,z(j)))/2; % Eddy diffusion @ j-1/2
    
    K(j,2)=(k(ustar,L,z(j+1))+k(ustar,L,z(j+2)))/2;% Eddy diffusion @ j+1/2
    
    alpha(j,1)=(2*K(j,1)*delx)/(u(ustar,z0,L,z(j+1))*(z(j+1)-z(j))*(z(j+2)-z(j)));
    
    alpha(j,2)=(2*K(j,2)*delx)/(u(ustar,z0,L,z(j+1))*(z(j+2)-z(j+1))*(z(j+2)-z(j)));

    % Constructing the Tridiagonal Matrix

    e(j)=-alpha(j,1);g(j)=-alpha(j,2);
    
    f(j)=1+alpha(j,1)+alpha(j,2);
    
    
end

vd=0.01;

beta=2*(z(2)-z(1))*vd/(k(ustar,L,z(1))+k(ustar,L,z(2)));


f(1)=f(1)-alpha(1,1)/(1+beta); % Boundary Condition at ground level

f(N-2)=f(N-2)-alpha(N-2,2); % Boundary Condition on surface layer height


% Source Boundary Condition 
N_s=find(z>=hs,1);

C(1,1:N)=0; % Source Condition

N_initial=1;

ubar=0;

for i=1:N_initial
    
    ubar=ubar+u(ustar,z0,L,z(N_s-i))+u(ustar,z0,L,z(N_s+i));
    
end

ubar=(ubar+u(ustar,z0,L,z(N_s)))/(2*N_initial+1);

C(1,N_s-N_initial:N_s+N_initial)=Q/ubar/(z(N_s+N_initial)-z(N_s-N_initial));

C(1,2)=C(1,1)*(1+beta);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Tridiagonal solver
for i=1:Nx-1
        
        b=C(i,2:N-1);
        
        C(i+1,2:N-1)=Tridiag_Solver(e,f,g,b);
        
        C(i+1,1)=C(i+1,2);
        
        C(i+1,N)=C(i+1,N-1);

        Z_bar(i)=sum(z(2:Nz).*C(i,:).*delz)/sum(C(i,:).*delz);
end

Cmass=C(5,:);
mass=0;

for j=1:Nz-2
    f1=Cmass(j)*u(ustar,z0,L,z(j));
    f2=Cmass(j+1)*u(ustar,z0,L,z(j+1));
    delz=z(j+1)-z(j);
    mass=mass+(f1+f2)*delz/2;
end


N_r=find(z>=1.5,1);

output=[C(51,N_r)/mass C(201,N_r)/mass C(801,N_r)/mass];

%  END OF PROGRAM

% Prairie Grass Run :
% Run 59 :  [C x z Z_bar]=surface_release(1000,0,@velocity,@eddy_diff,100,100,0.008,11,0.14);
% Run 26 :  [C x z Z_bar]=surface_release(1000,0,@velocity,@eddy_diff,100,100,0.008,-32,0.43); 
% i=101; j=8500 ; plot(C(i,2:j)/max(C(i,:)),z(2:j)/Z_bar(i)) ; ylim([0 4]);
% Grid System plot : X=[0:200:800]; for i=1:200; plot(X,z(i)*ones(1,length(X)),'k-'); hold on ;end


