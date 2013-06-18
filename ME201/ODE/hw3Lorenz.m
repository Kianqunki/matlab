function dx = hw3Lorenz(t,x)

global R_CONSTANT

dx = zeros(3,1); % a column vector
omega = 10;
b = 8./3;

dx(1) = omega*(x(2)-x(1));
dx(2) = R_CONSTANT*x(1) - x(2) - x(1)*x(3);
dx(3) = x(1)*x(2) - b*x(3);

end