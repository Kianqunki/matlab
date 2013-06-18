function [theta, C] = estimate_linear( X, U, params)

% dimension of observation and input
dim = 1;

% number of observations available
M = size(X,2);
p = params.p;
H = zeros(M, p);

% build matrix H
for i=1:M
    for j = min(i,p):-1:min(i,p)-p+1
        if ( j >= 1 )
            H(i,j) = U(dim,i-j+1);
        end
    end
end
% build the observations vector x as in the book
x = X';

% compute estimate of theta
temp = inv(H'*H);
theta = temp*H'*x;

% compute covariance matrix
C = params.sigma^2*temp;

% figure(1)
% spy(H)
% H
% keyboard;

end