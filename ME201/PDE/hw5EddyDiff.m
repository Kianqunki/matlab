function Kz = hw5EddyDiff(z, z0, L, u_star)
% compute the vertical eddy diffusivity Kz
zeta = z/L;

if ( L <= 0 )
    phi = (1-15*zeta)^(-0.25);
else
    phi = 1+4.7*zeta;
end
Km = 0.35*u_star*z/phi;

% alpha number
alpha = 1.35;
Kz = alpha*Km;

end