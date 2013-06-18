function hw2PrintFactor(Phi,NameList)
%hw1PrintFactor  Print the factor function on screen
%On input:
%Phi: The cell data structure that represents the factor function Phi

%List of measure values of the factor function
List=Phi.List;
%Number of values for each variable in the factor function
Val=Phi.Val;
%Name of the variables in the factor function (1 for X1, 2 for X2)
Variable=Phi.Variable;
%String array to convert numerical variable in Variable to text 
%(e.g. 1 -> 'a', 3 -> 'x3')
%VariableName = {'c' 'd' 't' 'g' 'i' 's' 'l' 'j'};
VariableName = NameList;
%VariableName = {'x1' 'x2' 'x3' 'x4'};

ass = zeros(1,length(Val));

for j=1:length(List)
    %Find the subscripts given the ValXYZ and index j
    %TODO: put this into a function hw1index2subs
    ndx = j;
    k = [1 cumprod(Val(1:end-1))];
    for i = length(Val):-1:1
        vi = rem(ndx-1, k(i)) + 1;
        vj = (ndx - vi)/k(i) + 1;
        ass(i) = vj;
        ndx = vi;
    end
    
    meas=hw2GetMeasure(Phi,ass);
    for i=1:length(Val)
        fprintf( 1, '%s=%d ', VariableName{Variable(i)}, ass(i)-1);
    end
    fprintf( 1, ': %f\n', meas );
end

fprintf( 1, '=================\n' );