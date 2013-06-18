A = hilb(50);
b = ones(50,1).*1/sqrt(50);
[U, S, V] = svd(A);
sDiag = diag(S);
S_inv = diag(1./sDiag);

% The original solution
x = V*S_inv*U'*b;

% The perturbed b vector
b(50) = b(50) - 1e-5;
% The new solution x'
xp = V*S_inv*U'*b;

% The relative change in the solution
change = norm(xp-x)/norm(x)

% Theoretical bound of the relative change
% Since ||dA|| = 0, ||db|| = 1e-5
bound = cond(A)*1e-5/norm(b)