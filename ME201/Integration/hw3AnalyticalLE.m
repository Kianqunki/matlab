function out = hw3Romberg(X,Y1,L,Xr,Yr,theta,q)
% initialization
Y2 = Y1+L;

% effective distance
xeff = Xr./cosd(theta);

%compute u1, u2
dem = sqrt(2).*sigmay(xeff);
u1 = ((Yr-Y1).*cosd(theta)-Xr.*sind(theta))./dem;
u2 = ((Yr-Y2).*cosd(theta)-Xr.*sind(theta))./dem;

% compute C
nom = q.*(erf(u1)-erf(u2));
dem = sqrt(2*pi).*cosd(theta).*sigmaz(xeff);

out = nom./dem;

end