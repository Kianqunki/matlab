function Phi = hw2FactorDivide(PhiXY, PhiY)
%hw2FactorDivide  Factor division

% %DEBUG
% NameList = {'c' 'd' 't' 'g' 'i' 's' 'l' 'j'};
% disp( 'Factor division' );
% hw2PrintFactor(PhiXY, NameList);
% hw2PrintFactor(PhiY,NameList);

% find the index of variable Y as in X
% a.k.a the "locations" of Y variable in X
[arr indexMarg dummy] = intersect(PhiXY.Variable,PhiY.Variable);

Phi.Variable = PhiXY.Variable;
Phi.Val = PhiXY.Val;

Phi.List = zeros(1,prod(Phi.Val));
SMALL = 0.00001;

for i=1:length(Phi.List)
    ass = hw2idx2sub(i,Phi.Val);
    assY = ass(indexMarg);
    idx2 = hw2sub2idx(assY, PhiY.Val);
    
    if ( abs(PhiXY.List(i)) < SMALL || abs(PhiY.List(idx2)) < SMALL)
        Phi.List(i) = 0;
    else
        Phi.List(i) = PhiXY.List(i)/PhiY.List(idx2);
    end
end
