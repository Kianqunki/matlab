function uout = hw5WindSpeed(z, z0, L, u_star)
% compute the horizontal wind speed
% paper 1 (by Nieuwstadt): equat. 2, 3
zeta = z/L;

if ( L <= 0 )
    x = (1-15*zeta)^0.25;
    psi = 2*log((1+x)/2)+log((1+x^2)/2)-2*atan(x)+pi/2;
else
    psi = -4.7*zeta;
end

uout = u_star/0.35*(log(z/z0)-psi);

end