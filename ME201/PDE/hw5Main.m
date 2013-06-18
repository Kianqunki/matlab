clear all;
close all;

% load data from Prairie Grass experiments
% Table 2 from the paper 2
load data;

% L and u* is given in the table2(:,2) and table2(:,3)
runNum = table2(:,1);
L = table2(:,2);
u_star = table2(:,3);
out50m = table2(:,4);
out200m = table2(:,5);
out800m = table2(:,6);
% parameters provided
% the surface roughness length z0
params.z0 = 0.008;
% height of the point source. paper 2, section 5
% h = 0.46;
params.h = 0.5;
% height of measurements. paper 2, section 5
params.zm = 1.5;
% maximum height of computation
params.zmax = 1500;
% maximum horizontal distance
params.xmax = 1000;

figure(2);
grid on
hold on
for k=1:length(runNum)
    [C, x, z] = hw5PDEsolver(@hw5WindSpeed, @hw5EddyDiff, L(k), u_star(k), params);
    zm = 1.5;
    [diff zIdx] = min(abs(zm-z));
    
    if ( ~isnan(out50m(k)) )
        xm = 50;
        [diff xIdx] = min(abs(xm-x));
%         keyboard;
        Cnew = hw5MassConv( C(xIdx,:), x, z, L(k), u_star(k), params )
%         C(xIdx,zIdx)
%         out50m(k)
        loglog(Cnew(zIdx),out50m(k), 'b*');
    end
    
    if ( ~isnan(out200m(k)) )
        xm = 200;
        [diff xIdx] = min(abs(xm-x));
        Cnew = hw5MassConv( C(xIdx,:), x, z, L(k), u_star(k), params );
%         C(xIdx,zIdx)
%         out200m(k)
        loglog(Cnew(zIdx),out200m(k), 'g*');
    end
    
    if ( ~isnan(out800m(k)) )
        xm = 800;
        [diff xIdx] = min(abs(xm-x));
        Cnew = hw5MassConv( C(xIdx,:), x, z, L(k), u_star(k), params );
%         C(xIdx,zIdx)
%         out800m(k)
        loglog(Cnew(zIdx),out800m(k), 'r*');
    end
end
V=axis;
X=[V(1) V(2)];
line(X,X);
line(X,2*X);
line(X,0.5*X);
