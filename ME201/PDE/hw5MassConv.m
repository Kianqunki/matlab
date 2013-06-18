function Cnew = hw5MassConv( C, x, z, L, u_star, params )
mass = 0;
Nz = length(z);
z0 = params.z0;

% compute all u and Kz
for j=1:Nz
    uVal(j) = hw5WindSpeed(z(j), z0, L, u_star);
end

for j=1:Nz-1
    f1 = C(j)*uVal(j);
    f2 = C(j+1)*uVal(j+1);
    delz = z(j+1)-z(j);
    mass = mass + (f1+f2)*0.5*delz;
end

Cnew = C/mass*1e3;

end