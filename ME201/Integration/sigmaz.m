function out = sigmaz(xr)
% sigma_z(x_r) in the equation

out = 0.14.*xr.*(1+0.0003*xr).^(-0.5);

end