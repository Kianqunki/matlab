function [Q, R] = QR_Gram(A)

m = size(A,1);
n = size(A,2);

Q = zeros(size(A));
R = zeros(n,n);

col = A(:,1);
Q(:,1) = col./norm(col);
R(1,1) = Q(:,1)'*col;
for i=2:n
    col = A(:,i);
    for j=1:i-1
        R(j,i) = Q(:,j)'*col;
        col = col - (Q(:,j)'*col)*Q(:,j);
    end
    Q(:,i) = col./norm(col);
    R(i,i) = Q(:,i)'*col;
end

end
