clear all;
load HW1.mat
[L, U] = my_lu(A,0);

% Ax = LUx = b => Ly = b; Ux = y
[nx, ny] = size(A);

% back-substitution to compute y
y = zeros(nx, 1);
y(1) = b(1)/L(1,1);
for k=2:nx
    sum = L(k,1:k-1)*y(1:k-1);
    y(k) = (b(k)-sum)/L(k,k);
end
y

% back-substitution to compute x
x = zeros(nx,1);
x(nx) = y(nx)/U(nx,nx);
for k=nx-1:-1:1
    sum = U(k,k+1:nx)*x(k+1:nx);
    x(k) = (y(k)-sum)/U(k,k);
end
x

clear all;
load HW1.mat
[L, U] = my_lu(A,0);

% Ax = LUx = b => Ly = b; Ux = y
[nx, ny] = size(A);

% back-substitution to compute y
y = zeros(nx, 1);
y(1) = b(1)/L(1,1);
for k=2:nx
    sum = L(k,1:k-1)*y(1:k-1);
    y(k) = (b(k)-sum)/L(k,k);
end
y

% back-substitution to compute x
x = zeros(nx,1);
x(nx) = y(nx)/U(nx,nx);
for k=nx-1:-1:1
    sum = U(k,k+1:nx)*x(k+1:nx);
    x(k) = (y(k)-sum)/U(k,k);
end
x