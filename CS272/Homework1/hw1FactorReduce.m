function PhiRed = hw1FactorReduce(Phi, Vars_U, Value_u)
%hw1FactorReduce  Print the reduced factor function on screen
%On input:
%Phi: The cell data structure that represents the factor function Phi
%Vars_U: The vector that stores the numerical name of variables of
%context U
%Value_u: The vector that stores the values u of variables of context U
%On output:
%PhiRed: The cell data structure that represents the reduced factor
%function

%List of measure values of the factor function
List=Phi{1};
%Number of values for each variable in the factor function
Val=Phi{2};
%Name of the variables in the factor function (1 for X1, 2 for X2)
Variable=Phi{3};

% the index of variables in U as in X
indexU=zeros(1,length(Vars_U));
% the index of other variables in X
indexRed=zeros(1,length(Val)-length(Vars_U)); 
ndx=1;
idx=1;

% find the index of variables in U as in X
% a.k.a the "locations" of U variables in X
for i=1:length(Variable)
    index = find(Vars_U==Variable(i));
    if (size(index)>0)
        indexU(ndx) = i;
        ndx=ndx+1;  %not the best method
    else
        indexRed(idx)=i;
        idx=idx+1;
    end
end

ass = zeros(1,length(indexRed));
%initialize the reduced factor function
ValRed = Val(indexRed);
VariableRed = Variable(indexRed);
ListRed = zeros(1,prod(ValRed));

% pre-compute the index due to fixed values in u
idx = 1;
k1 = [1 cumprod(Val(1:end-1))];
for i = 1:length(Vars_U),
    idx = idx + (Value_u(i)-1)*k1(indexU(i));
end

% find all the measure values for the reduced factor function
for j=1:length(ListRed)
    %Find the subscripts given the ValRed and index j
    %Not the best way to generate subscript
    %TODO: put this into a function hw1index2subs
    ndx = j;
    k2 = [1 cumprod(ValRed(1:end-1))];
    for i = length(ValRed):-1:1
        vi = rem(ndx-1, k2(i)) + 1;
        vj = (ndx - vi)/k2(i) + 1;
        ass(i) = vj;
        ndx = vi;
    end
    
    %Find the index in the original factor function
    idx2 = idx;
    for i = 1:length(indexRed)
        idx2 = idx2 + (ass(i)-1)*k1(indexRed(i));
    end
    
    ListRed(j) = List(idx2);
end

PhiRed = {ListRed ValRed VariableRed};