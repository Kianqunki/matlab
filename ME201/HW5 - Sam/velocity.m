function u=velocity(ustar,z0,L,z)

if L>=0
    Psi=-4.7*z/L;
else
    x=(1-15*z/L)^0.25;
    Psi=2*log((1+x)/2)+log((1+x^2)/2)-2*atan(x)+pi/2;
end


u=(ustar/0.4)*[log(z/z0)-Psi];

%End OF PROGRAM

