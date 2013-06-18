function out = sigmay(xr)
% sigma_y(x_r) in the equation

out = 0.16.*xr.*(1+0.0004*xr).^(-0.5);

end