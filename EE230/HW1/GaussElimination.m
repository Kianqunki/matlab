function x = GaussElimination(A,b)
% % solving linear equation Ax=b by Gaussian elimination
% % Example:
% A=[1 2 4; 5 8 9; 3 7 5];
% b = [2 8 9]';
% GaussElimination(A,b)

[nx, ny] = size(A);
aug = [A,b];
x=zeros(size(b));

for k=1:nx-1
    [cm, im] = max(abs(aug(k:nx,k)));
    if cm ~= 0
        % swap rows
        dum = aug(k,:);
        aug(k,:) = aug(im+k-1,:);
        aug(im+k-1,:) = dum;
        
        for i=k+1:nx
            mult = aug(i,k)/aug(k,k);
            aug(i,k:nx+1) = aug(i,k:nx+1) - mult.*aug(k,k:nx+1);
        end
    end
end

AA = aug(:,1:nx);
bb = aug(:,nx+1);

x(nx) = bb(nx)/AA(nx,nx);
for k=nx-1:-1:1
    sum = AA(k,k+1:nx)*x(k+1:nx);
    x(k) = (bb(k)-sum)/AA(k,k);
end

end