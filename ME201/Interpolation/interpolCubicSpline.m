function yout = interpolCubicSpline(x, y, xin)
% Input: 
% x, y: table of values. N-by-1 vectors.
% xin: the value at which the function is interpolated
% Output:
% yout: the interpolated value

% Optional error checking: 
% 1) x, y has the same length? More than some minimum number?
% 2) No identical values, esp. in x.
% assuming natural spline

% find the index of value in the table closest to xin
[diffVal, idx1] = min(abs(x-xin));
if (x(idx1) < xin)
    idx2 = idx1+1;
else
    idx2 = idx1;
    idx1 = idx1-1;
end
%x(idx1), x(idx2) are the two closest values to xin

N = length(x);

% Compute A, B, C, D
A = (x(idx2)-xin)./(x(idx2)-x(idx1));
B = 1-A;
factor = (x(idx2)-x(idx1)).^2/6;
C = factor*(A^3-A);
D = factor*(B^3-B);

% Compute e, f, g, b
% initialize
e = zeros(N,1);
f = zeros(N,1);
g = zeros(N,1);
b = zeros(N,1);

diff_x = diff(x)./6;
diff_y = diff(y)./6;
diff_xy = diff_y./diff_x;
% diff(X) returns [X(2)-X(1) X(3)-X(2) ... X(n)-X(n-1)]
b(2:N-1) = diff(diff_xy);
e(2:N-1) = diff_x(1:end-1);
g(2:N-1) = diff_x(2:end);
f(2:N-1) = 2*(diff_x(1:end-1)+diff_x(2:end));
% Boundary condition for natural spline
f(1) = 1;
f(N) = 1;
% For natural splines, other entries in b, e, g are zeros

% Compute y(2)
y2 = Tridiag_Solver(e,f,g,b);

% Compute yout
yout = A*y(idx1)+B*y(idx2)+C*y2(idx1)+D*y2(idx2);