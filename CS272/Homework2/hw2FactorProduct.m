function Product = hw2FactorProduct( PhiXY, PhiYZ )
%hw1FactorProduct  Return the factor product of factor PhiXY and factor YZ,
%given the common part Y.
%On input:
%PhiXY: The cell data structure that represents the factor function Phi1
%PhiYZ: The cell data structure that represents the factor function Phi2
%On output:
%Product: The cell data structure that represents the product factor function

%Difference from hw1FactorProduct:
%No user-specified X,Y,Z (especially shared variables Y)

%List of measure values of the factor function
ListXY=PhiXY.List;
%Number of values for each variable in the factor function
ValXY=PhiXY.Val;     
%Name of the variables in the factor function (1 for X1, 2 for X2)
VariableXY = PhiXY.Variable;    

%List of measure values of the factor function
ListYZ=PhiYZ.List;
%Number of values for each variable in the factor function
ValYZ=PhiYZ.Val;
%Name of the variables in the factor function (1 for X1, 2 for X2)
VariableYZ = PhiYZ.Variable;     

%Name of the variables in the factor product
VariableXYZ = unique([VariableXY VariableYZ]);
ValXYZ = zeros(1, length(VariableXYZ) );

mapXY = zeros(1,length(VariableXY)); %index of variables in XY in factor XYZ
mapYZ = zeros(1,length(VariableYZ)); %index of variables in YZ in factor XYZ

% find the indexes of variables of X, Y in Z
% and fill the ValXYZ
for i=1:length(ValXYZ)
    index1 = find(VariableXY==VariableXYZ(i));
    index2 = find(VariableYZ==VariableXYZ(i));
    
    if (isempty(index1))
        % not in PhiXY, aka Z
        ValXYZ(i) = ValYZ(index2);
        mapYZ(index2) = i;
    elseif (isempty(index2))
        % not in PhiYZ, aka X
        ValXYZ(i) = ValXY(index1);
        mapXY(index1) = i;
    elseif (VariableXY(index1) == VariableYZ(index2))
        % that means Y
        ValXYZ(i) = ValXY(index1);
        mapXY(index1) = i;
        mapYZ(index2) = i;
    else
        % ERROR: error handling here
    end
end

% initialize components of the product factor function
ListXYZ = zeros(1,prod(ValXYZ));

for j=1:length(ListXYZ)
    %Find the subscripts given the ValXYZ and index j
    assXYZ = hw2idx2sub(j, ValXYZ);
    %Rearrange the subscripts to find subsripts for PhiXY and PhiYZ
    assXY = assXYZ(mapXY);
    assYZ = assXYZ(mapYZ);

    meas1 = hw2GetMeasure(PhiXY,assXY);
    meas2 = hw2GetMeasure(PhiYZ,assYZ);
    ListXYZ(j) = meas1*meas2;
end

Product.List = ListXYZ;
Product.Val = ValXYZ;
Product.Variable = VariableXYZ;
