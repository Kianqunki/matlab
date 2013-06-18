function [Q, R] = QR_Householder(A)

m = size(A,1);
n = size(A,2);

Q = eye(m,m);
R = A;

for i=1:n
    col = R(i:end,i);
    e1 = zeros(size(col));
    e1(1) = 1;
    v = col - norm(col)*e1;
    subH = eye(length(col)) - 2*v*v'/(v'*v);
    
    Qi = eye(m,m);
    Qi(i:end,i:end) = subH;
    R = Qi*R;
    Q = Q*Qi';
end

Q = Q(:,1:n);
R = R(1:n,:);
end
