function PhiMarg = hw2FactorMarginalize(Phi, Var_Y)
%hw1FactorReduce  Normalize a factor, given a variable
%On input:
%Phi: The cell data structure that represents the factor function Phi.
%Var_Y: The vector that stores the numerical name of variable Y.
%PhiMarg: The cell data structure that represents the marginalized factor
%function

%Difference with hw1FactorMarginalize:
%More general Var_Y: empty, scalar, vector. To use in clique tree.

%List of measure values of the factor function
List=Phi.List;
%Number of values for each variable in the factor function
Val=Phi.Val;
%Name of the variables in the factor function (1 for X1, 2 for X2)
Variable=Phi.Variable;

%check for empty Var_Y
if ( isempty(Variable) || isempty(Var_Y) )
    PhiMarg = Phi;
    return;
end

% find the index of variable Y as in X
% a.k.a the "locations" of Y variable in X
[VariableMarg indexMarg] = setdiff(Variable, Var_Y);

%initialize the marginalized factor function
ValMarg = Val(indexMarg);
ListMarg = zeros(1,prod(ValMarg));


% find all the measure values for the reduced factor function
for j=1:length(List)
    %Find the subscripts given the ValMarg and index j
    ass = hw2idx2sub(j,Val);
    assMarg = ass(indexMarg);
    idx2 = hw2sub2idx(assMarg, ValMarg);
    
    ListMarg(idx2) = ListMarg(idx2) + List(j);
end

PhiMarg.List = ListMarg;
PhiMarg.Val = ValMarg;
PhiMarg.Variable = VariableMarg;