clear all
close all

%define factor function Phi1
%List of measure values of the factor function
ListAB = [30 1 5 10];
%Number of values for each variable in the factor function
ValAB = [2 2];         
%Name of the variables in the factor function (1 for X1, 2 for X2)
VariableAB = [1 2];    
PhiAB = {ListAB ValAB VariableAB};
hw1PrintFactor(PhiAB);

%define factor function Phi2
ListBC = [100 1 1 100]; 
ValBC = [2 2];        
VariableBC = [2 3]; 
PhiBC = {ListBC ValBC VariableBC};
hw1PrintFactor(PhiBC);

%define factor function Phi3
ListCD = [1 100 100 1];
ValCD = [2 2];        
VariableCD = [3 4]; 
PhiCD = {ListCD ValCD VariableCD};
hw1PrintFactor(PhiCD);

%define factor function Phi4
ListDA = [100 1 1 100];
ValDA = [2 2];    
VariableDA = [4 1];
PhiDA = {ListDA ValDA VariableDA};
hw1PrintFactor(PhiDA);

%find the factor product
PhiABC = hw1FactorProduct(PhiAB, PhiBC, [1],[2],[3]);
PhiCDA = hw1FactorProduct(PhiCD, PhiDA, [3],[4],[1]);
PhiABCD = hw1FactorProduct(PhiABC, PhiCDA, [2],[1 3],[4]);
hw1PrintFactor(PhiABCD);

%Create the normalized factor product
disp( 'Check factor normalization' );
NormListABCD = PhiABCD{1}./sum(PhiABCD{1});
NormPhiABCD = {NormListABCD PhiABCD{2} PhiABCD{3}};
hw1PrintFactor(NormPhiABCD);

disp( 'Check factor marginalization' );
disp( 'Results in page 105' );
NormPhiBCD = hw1FactorMarginalize(NormPhiABCD, [1]);
NormPhiBD = hw1FactorMarginalize(NormPhiBCD, [3]);
NormPhiB = hw1FactorMarginalize(NormPhiBD, [4]);
hw1PrintFactor(NormPhiB);

% for Example 4.2
disp( 'Results in Example 4.2' );
NormPhiABC = hw1FactorMarginalize(NormPhiABCD, [4]);
NormPhiAB = hw1FactorMarginalize(NormPhiABC, [3]);
hw1PrintFactor(NormPhiAB);
