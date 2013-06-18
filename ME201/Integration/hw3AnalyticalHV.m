function out = hw3AnalyticalHV(X,Y1,L,Xr,Yr,thetad,q)
% initialization
Y2 = Y1+L;
theta = thetad;

% effective distance
xeff = Xr./cosd(theta);
x1 = Xr.*cosd(theta) + (Yr-Y1).*sind(theta);
x2 = Xr.*cosd(theta) + (Yr-Y2).*sind(theta);

%compute u1, u2
dem1 = sqrt(2).*sigmay(x1);
t1 = ((Yr-Y1).*cosd(theta)-Xr.*sind(theta))./dem1;
if (x2 <= 0)
    t2 = -sign(theta)*inf;
else
    dem2 = sqrt(2).*sigmay(x2);
    t2 = ((Yr-Y2).*cosd(theta)-Xr.*sind(theta))./dem2;
end

% compute C
nom = q.*(erf(t1)-erf(t2));
dem = sqrt(2*pi).*cosd(theta).*sigmaz(xeff);

out = nom./dem;

end