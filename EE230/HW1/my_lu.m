function [L, U] = my_lu(A, flag)

[nx, ny] = size(A);

if (flag == 1)
    % U and unit L
    % U is a matrix of size of A
    % L is a square matrix nx by nx
    L = eye(nx, nx);
    
    for k=1:nx-1
        % the inverse of product of Eij's
        Ln = eye(nx,nx);
        
        for i=k+1:nx
            % multiplier
            l = -A(i,k)/A(k,k);
            % subtract the i'th row by the k'th row multiplied by l
            A(i,k:ny) = A(i,k:ny) + l*A(k,k:ny);
            % construct Eij inversed
            Ln(i,k) = -l;
        end
        
        L = L*Ln;
    end
    
    U = A;
elseif (flag == 0)
    % Unit U and L
    
    % U is a matrix of size of A
    % L is a square matrix nx by nx
    L = eye(nx, nx);
    
    for k=1:nx
        % the inverse of product of Eij's
        Ln = eye(nx,nx);
        % the inverse of Ei
        Lm = eye(nx,nx);
        % construct Ei inversed
        Lm(k,k) = A(k,k);
        % divided the row by the pivot
        A(k,:) = A(k,:)/A(k,k);
        
        for i=k+1:nx           
            % construct Eij inversed.
            Ln(i,k) = A(i,k);
            % subtract the i'th row by the k'th row multiplied by l
            A(i,k:ny) = A(i,k:ny) - A(i,k)*A(k,k:ny);
        end
        
        L = L*Lm*Ln;
    end
    
    U = A;
else
    % Wrong input. Do  nothing. Warning if needed.
end

end