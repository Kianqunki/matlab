function out = hw3Gaussian(X,Y1,L,Xr,Yr,theta,q,order)
% initialization
Y2 = Y1+L;

% Compute Legendre polynomials iteratively
p(1,1)  =1;     % P0 = 1
p(2,1:2)=[1 0]; % P1 = x
for k=2:order
    % entry k is for P(k-1)
    p(k+1,1:k+1)=((2*k-1)*[p(k,1:k) 0]-(k-1)*[0 0 p(k-1,1:k-1)])/k;
end

Pn = p(order+1,:);
Pn_deriv = polyder(Pn);

% find the zeros for the Pn
x = roots(Pn);
% find the weights: from the formula provided on Wikipedia
w = 2./((1-x.^2).*(polyval(Pn_deriv,x)).^2);

% Shift from [-1,1] to [Y1,Y2]:
Y = (Y2+Y1)/2 + L.*x/2;
w = L/2 * w;

% evaluate the function at each point
f = hw3SourceContrib(X,Y,Xr,Yr,theta,q);
% dot products of w and f
out = dot(w,f);
end