% script to test GuassElimination function

for i=1:10
    A = rand(i);
    b = rand(i,1);
    X = GaussElimination(A,b)
    X2 = Gauss_Elimination(A,b)
    
    % compare with the standard solution
    Xtrue = A\b
    input('Wait');
end