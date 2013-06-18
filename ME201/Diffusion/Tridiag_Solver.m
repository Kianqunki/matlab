function X = Tridiag_Solver(e, f, g, b)
N = length(f);

X = zeros(size(b));
%Gauss elimination
for k = 2:N
    mult = e(k)/f(k-1);
    f(k) = f(k) - mult*g(k-1);
    b(k) = b(k) - mult*b(k-1);
end

X(N) = b(N)/f(N);
for k = N-1:-1:1
    X(k) = (b(k) - g(k)*X(k+1))/f(k);
end

end