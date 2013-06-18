function Product = hw1FactorProduct( PhiXY, PhiYZ, X, Y, Z)
%hw1FactorProduct  Return the factor product of factor PhiXY and factor YZ,
%given the common part Y.
%On input:
%PhiXY: The cell data structure that represents the factor function Phi1
%PhiYZ: The cell data structure that represents the factor function Phi2
%X, Y, Z: The vectors that specify the variables that are in X, Y, Z
%On output:
%Product: The cell data structure that represents the product factor function

%List of measure values of the factor function
ListXY=PhiXY{1};
%Number of values for each variable in the factor function
ValXY=PhiXY{2};     
%Name of the variables in the factor function (1 for X1, 2 for X2)
VariableXY = PhiXY{3};    

%List of measure values of the factor function
ListYZ=PhiYZ{1};
%Number of values for each variable in the factor function
ValYZ=PhiYZ{2};
%Name of the variables in the factor function (1 for X1, 2 for X2)
VariableYZ = PhiYZ{3};     

%TODO: check if XYZ is compatiable with PhiXY, PhiYZ.
%check sum(VariableXY) == sum(X) + sum(Y)
%check sum(VariableYZ) == sum(Y) + sum(Z)

indexX = zeros(length(X),1);        %index of variables in X in factor XY
indexY_XY = zeros(length(Y),1);     %index of variables in Y in factor XY
indexY_YZ = zeros(length(Y),1);     %index of variables in Y in factor YZ
indexZ = zeros(length(Z),1);        %index of variables in Z in factor YZ

%find the indexes
for i=1:length(X)
    indexX(i) = find(VariableXY==X(i));
end

for i=1:length(Y)
    indexY_XY(i) = find(VariableXY==Y(i));
    indexY_YZ(i) = find(VariableYZ==Y(i));
end

for i=1:length(Z)
    indexZ(i) = find(VariableYZ==Z(i));
end

ValX = ValXY(indexX);
ValY = ValXY(indexY_XY);
ValZ = ValYZ(indexZ);

% initialize components of the product factor function
ListXYZ = zeros(1,prod(ValX)*prod(ValY)*prod(ValZ));
ValXYZ = [ValX ValY ValZ];
VariableXYZ = [X Y Z];

assXYZ = zeros(1,length(VariableXYZ));
assXY = zeros(1, length(VariableXY));
assYZ = zeros(1, length(VariableYZ));

for j=1:length(ListXYZ)
    %Find the subscripts given the ValXYZ and index j
    ndx = j;
    k = [1 cumprod(ValXYZ(1:end-1))];
    for i = length(ValXYZ):-1:1
        vi = rem(ndx-1, k(i)) + 1;
        vj = (ndx - vi)/k(i) + 1;
        assXYZ(i) = vj;
        ndx = vi;
    end
    
    %Rearrange the subscripts to find subsripts for PhiXY and PhiYZ
    for i=1:length(X)
        assXY(indexX(i))=assXYZ(i);
    end
    for i=1:length(Y)
        assXY(indexY_XY(i))=assXYZ(length(X)+i);
        assYZ(indexY_YZ(i))=assXYZ(length(X)+i);
    end
    for i=1:length(X)
        assYZ(indexZ(i))=assXYZ(length(VariableXY)+i);
    end
    
    meas1 = hw1GetMeasure(PhiXY,assXY);
    meas2 = hw1GetMeasure(PhiYZ,assYZ);
    ListXYZ(j) = meas1*meas2;    
end

Product = {ListXYZ ValXYZ VariableXYZ};
