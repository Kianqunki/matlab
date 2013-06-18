function out = hw3Romberg(X,Y1,L,Xr,Yr,theta,q,order)
% initialization
Y2 = Y1+L;
I  = zeros(order,order);

% compute the base integrals by trapezoidal method
for i=1:order
    N = 2^i;
    dY = L./N;
    Y = Y1:dY:Y2;
    % thus length(Y) should be N+1;
    
%     f = zeros(length(Y),1);    
%     for j = 1:length(Y)
%         f(j) = hw3SourceContrib(X,Y(j),Xr,Yr,dY,theta,q);
%         I(i) = I(i) + f(j);
%     end
    
    f = hw3SourceContrib(X,Y,Xr,Yr,theta,q);
    I(i,1) = I(i) + sum(f(2:N));
    I(i,1) = I(i) + 0.5*(f(1)+f(N+1));
    I(i,1) = I(i)*dY;
end

% compute the final integral by Romberg method
for j=1:(order-1)
    num = 2^j;
%     keyboard;
    for i=1:(order-j)
        I(i,j+1) = num.*I(i+1,j) - I(i,j);
        I(i,j+1) = I(i,j+1)./(num-1);
    end
end

out = I(1,order);

% compute Romberg by using trapezoidal and interpolation.