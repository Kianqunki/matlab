clear all
close all

%Purpose: check factor marginalization function

%define factor function Phi1
%List of measure values of the factor function
ListAB = [0.5 0.1 0.3 0.8 0 0.9];  
%Number of values for each variable in the factor function
ValAB = [3 2];       
%Name of the variables in the factor function (1 for X1, 2 for X2)
VariableAB = [1 2];
PhiAB = {ListAB ValAB VariableAB};
hw1PrintFactor(PhiAB);

%define factor function Phi2
ListBC = [0.5 0.1 0.7 0.2]; 
ValBC = [2 2];      
VariableBC = [2 3];  
PhiBC = {ListBC ValBC VariableBC};
hw1PrintFactor(PhiBC);

%find the factor product
PhiABC = hw1FactorProduct(PhiAB, PhiBC, [1],[2],[3]);
hw1PrintFactor(PhiABC);
PhiAC = hw1FactorMarginalize(PhiABC, [2] );
hw1PrintFactor(PhiAC);