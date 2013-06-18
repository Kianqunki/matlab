function PhiMarg = hw1FactorReduce(Phi, Var_Y)
%hw1FactorReduce  Print the reduced factor function on screen
%On input:
%Phi: The cell data structure that represents the factor function Phi.
%Var_Y: The vector that stores the numerical name of variable Y.
%PhiMarg: The cell data structure that represents the marginalized factor
%function

%List of measure values of the factor function
List=Phi{1};
%Number of values for each variable in the factor function
Val=Phi{2};
%Name of the variables in the factor function (1 for X1, 2 for X2)
Variable=Phi{3};

% find the index of variable Y as in X
% a.k.a the "locations" of Y variable in X

% the index of variable Y as in X
indexY = find(Variable==Var_Y); 
% the index of other variables in X
indexMarg = [1:(indexY-1) (indexY+1):length(Variable)]; 

ass = zeros(1,length(indexMarg));
%initialize the marginalized factor function
ValMarg = Val(indexMarg);
VariableMarg = Variable(indexMarg);
ListMarg = zeros(1,prod(ValMarg));

k1 = [1 cumprod(Val(1:end-1))];

% find all the measure values for the reduced factor function
for j=1:length(ListMarg)
    %Find the subscripts given the ValMarg and index j
    %Not the best way to generate subscript
    %TODO: put this into a function hw1index2subs
    ndx = j;
    k2 = [1 cumprod(ValMarg(1:end-1))];
    for i = length(ValMarg):-1:1
        vi = rem(ndx-1, k2(i)) + 1;
        vj = (ndx - vi)/k2(i) + 1;
        ass(i) = vj;
        ndx = vi;
    end
    
    %Find the index in the original factor function
    idx = 1;
    for i = 1:length(indexMarg)
        idx = idx + (ass(i)-1)*k1(indexMarg(i));
    end
    
    %Summing over all values of Y
    for i = 1:Val(indexY)
        %Find the index for the given value
        idx2 = idx + (i-1)*k1(indexY);
        %Retrieve the values and sum up
        ListMarg(j) = ListMarg(j) + List(idx2); 
    end
end

PhiMarg = {ListMarg ValMarg VariableMarg};