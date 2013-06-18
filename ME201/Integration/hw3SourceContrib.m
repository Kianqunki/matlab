function contrib = hw3SourceContrib(X,Y,Xr,Yr,theta,q)
% for simplicity and clarity, the function is SISO and in long form

% co-ordinate transformation
xr = Xr.*cosd(theta) + (Yr-Y).*sind(theta);
yr = (Yr-Y).*cosd(theta) - Xr.*sind(theta);

% compute sigma functions' outputs
% % sigma_y(x_r) in the equation
% sigma2 = 0.16.*xr.*(1+0.0004*xr).^(-0.5);
% % sigma_y^2(x_r) in the equation
% %sigma1 = (xr.*2).^2;
% sigma1 = sigma2.^2;
% % sigma_z(x_r) in the equation
% sigma3 = 0.14.*xr.*(1+0.0003*xr).^(-0.5);

% compute sigma functions' outputs
sigma2 = sigmay(xr);
sigma1 = sigma2.^2;
sigma3 = sigmaz(xr);

dem = pi.*sigma2.*sigma3;

contrib = exp(-yr.^2./(2.*sigma1));
contrib = (q.*contrib)./dem;