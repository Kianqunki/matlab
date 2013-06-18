clear all;
close all;

NameList = {'w' 'l' 'a' 't' 'r' 'c' 's' 'p'};

% use struct to define factors instead of cell array
% as later, when an array of factors is defined, cell array reference is
% confusing

%Weather
%List of measure values of the factor function
PhiW.List = [0.2 0.8];
%Number of values for each variable in the factor function
PhiW.Val = [2];         
%Name of the variables in the factor function (1 for X1, 2 for X2)
PhiW.Variable = [1];    
%PhiC = {ListC ValC VariableC};
%hw1PrintFactor(PhiC,NameList);

%Location
PhiL.List = [0.6 0.4];
PhiL.Val = [2];
PhiL.Variable = [2];    
%PhiD = {ListD ValD VariableD};
%hw1PrintFactor(PhiD,NameList);

%Advertising
PhiA.List = [0.8 0.2];
PhiA.Val = [2];
PhiA.Variable = [3];
%PhiI = {ListI ValI VariableI};
%hw1PrintFactor(PhiI,NameList);

%Time
PhiT.List = [0.3 0.7];
PhiT.Val = [2];
PhiT.Variable = [4];
%PhiS = {ListS ValS VariableS};
%hw1PrintFactor(PhiS,NameList);

%Traffic
PhiR.List = [0.9 0.6 0.5 0.2 0.1 0.4 0.5 0.8];
PhiR.Val = [2 2 2];
PhiR.Variable = [1 2 5];
%PhiT = {ListT ValT VariableT};
%hw2PrintFactor(PhiR,NameList);

%Cost
PhiC.List = [0.9 0.4 0.7 0.2 0.1 0.6 0.3 0.8];
PhiC.Val = [2 2 2];
PhiC.Variable = [3 4 6];
%PhiL = {ListL ValL VariableL};
%hw2PrintFactor(PhiC,NameList);

%Sale
PhiS.List = [0.9 0.8 0.75 0.7 0.3 0.4 0.2 0.1 0.1 0.2 0.25 0.3 0.7 0.6 0.8 0.9];
PhiS.Val = [2 2 2 2];
PhiS.Variable = [3 4 5 7];
%PhiG = {ListG ValG VariableG};
% hw2PrintFactor(PhiS,NameList);

%Profit
PhiP.List = [0.1 0.7 0.2 0.2 0.8 0.2 0.3 0.6 0.1 0.1 0.5 0.2];
PhiP.Val = [2 2 3];
PhiP.Variable = [6 7 8];
%PhiJ = {ListJ ValJ VariableJ};
%hw2PrintFactor(PhiP,NameList);

%Set of factors
Factors = [PhiW PhiL PhiA PhiT PhiR PhiC PhiS PhiP];

%Define clique tree
Tree.varName = [1 2 3 4 5 6 7 8];
Tree.Val = [2 2 2 2 2 2 2 3];
Tree.nodes = {[1 2 5] [3 4 5 6 7] [6 7 8]};
% represent edges by using NxN matrix. Any better idea?
% count the edge only one
Tree.edges = [0 1 0;
    0 0 1;
    0 0 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Testing zone

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calibrate the tree and compute the marginals for each variable
[Beta, Mu] = hw2CalibrateCliqueTree(Tree, Factors);
disp( 'Printout all beliefs Beta_i' );
for i=1:length(Beta)
    hw2PrintFactor(Beta(i),NameList);
end

% %Calculate the marginals
% P_1 = hw2FactorMarginalize( Beta(1), [2 5] );
% P_2 = hw2FactorMarginalize( Beta(1), [1 5] );
% P_3 = hw2FactorMarginalize( Beta(2), [4 5 6 7] );
% P_4 = hw2FactorMarginalize( Beta(2), [3 5 6 7] );
% P_5 = hw2FactorMarginalize( Beta(2), [3 4 6 7] );
% P_6 = hw2FactorMarginalize( Beta(2), [3 4 5 7] );
% P_7 = hw2FactorMarginalize( Beta(3), [6 8] );
% P_8 = hw2FactorMarginalize( Beta(3), [6 7] );
% 
% disp( 'Printout all marginals for variables' );
% hw2PrintFactorNormalized(P_1,NameList);
% hw2PrintFactorNormalized(P_2,NameList);
% hw2PrintFactorNormalized(P_3,NameList);
% hw2PrintFactorNormalized(P_4,NameList);
% hw2PrintFactorNormalized(P_5,NameList);
% hw2PrintFactorNormalized(P_6,NameList);
% hw2PrintFactorNormalized(P_7,NameList);
% hw2PrintFactorNormalized(P_8,NameList);

%Automatically
% for each variable in tree
disp( 'Printout all marginals for variables' );
P = repmat(struct('List', [], 'Val', [], 'Variable', []), length(Tree.varName), 1);
for i=1:length(Tree.varName)
    idx=0;
    for j=1:length(Beta)
        if (any(Beta(j).Variable == Tree.varName(i) ) )
            idx=j;
            break;
        end
    end
    
    P(i) = hw2FactorMarginalize( Beta(idx), setdiff(Beta(idx).Variable,Tree.varName(i)) );
    hw2PrintFactor(P(i),NameList);
end
