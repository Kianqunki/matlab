clear all;
close all;

NameList = {'c' 'd' 't' 'g' 'i' 's' 'l' 'j'};

% use struct to define factors instead of cell array
% as later, when an array of factors is defined, cell array reference is
% confusing

%PhiC
%List of measure values of the factor function
PhiC.List = [0.5 0.5];
%Number of values for each variable in the factor function
PhiC.Val = [2];         
%Name of the variables in the factor function (1 for X1, 2 for X2)
PhiC.Variable = [1];    
%PhiC = {ListC ValC VariableC};
%hw1PrintFactor(PhiC,NameList);

%PhiD.
PhiD.List = [0.4 0.8 0.6 0.2];
PhiD.Val = [2 2];
PhiD.Variable = [1 2];    
%PhiD = {ListD ValD VariableD};
%hw1PrintFactor(PhiD,NameList);

%PhiI.
PhiI.List = [0.6 0.4];
PhiI.Val = [2];
PhiI.Variable = [5];
%PhiI = {ListI ValI VariableI};
%hw1PrintFactor(PhiI,NameList);

%PhiS.
PhiS.List = [0.95 0.2 0.05 0.8];
PhiS.Val = [2 2];
PhiS.Variable = [3 6];
%PhiS = {ListS ValS VariableS};
%hw1PrintFactor(PhiS,NameList);

%PhiT.
PhiT.List = [0.9 0.4 0.1 0.6];
PhiT.Val = [2 2];
PhiT.Variable = [5 3];
%PhiT = {ListT ValT VariableT};
%hw1PrintFactor(PhiT,NameList);

%PhiL.
PhiL.List = [0.1 0.4 0.99 0.9 0.6 0.01];
PhiL.Val = [3 2];
PhiL.Variable = [4 7];
%PhiL = {ListL ValL VariableL};
%hw1PrintFactor(PhiL,NameList);

%PhiG.
PhiG.List = [0.3 0.05 0.9 0.5 0.4 0.25 0.08 0.3 0.3 0.7 0.02 0.2];
PhiG.Val = [2 2 3];
PhiG.Variable = [2 3 4];
%PhiG = {ListG ValG VariableG};
%hw1PrintFactor(PhiG,NameList);

%PhiJ.
PhiJ.List = [0.9 0.4 0.3 0.1 0.1 0.6 0.7 0.9];
PhiJ.Val = [2 2 2];
PhiJ.Variable = [6 7 8];
%PhiJ = {ListJ ValJ VariableJ};
%hw1PrintFactor(PhiG,NameList);

%Set of factors
Factors = [PhiC PhiD PhiT PhiG PhiI PhiS PhiL PhiJ];

%Define clique tree
Tree.varName = [1 2 3 4 5 6 7 8];
Tree.Val = [2 2 2 3 2 2 2 2];
Tree.nodes = {[1 2] [4 2 3] [4 3 6] [5 3] [6 7 4] [8 7 6]};
% represent edges by using NxN matrix. Any better idea?
% Tree.edges = [0 1 0 0 0 0;
%     1 0 1 0 0 0;
%     0 1 0 1 1 0;
%     0 0 1 0 0 0;
%     0 0 1 0 0 1;
%     0 0 0 0 1 0];
% count the edge only one
Tree.edges = [0 1 0 0 0 0;
    0 0 1 0 0 0;
    0 0 0 1 1 0;
    0 0 0 0 0 0;
    0 0 0 0 0 1;
    0 0 0 0 0 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Testing zone
% length(Tree.nodes);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% idx = find( Tree.edges == 1);
% [i j] = ind2sub([N N], idx);
% Mu = repmat(struct('List', [], 'Val', [], 'Variable', []), N, N);
% for k=1:length(idx)
%     [arr, idx1, idx2] = intersect(Tree.nodes{i(k)},Tree.nodes{j(k)});
%     Mu(idx(k)).Variable = arr;
%     Mu(idx(k)).Val = Tree.Val(arr);
%     Mu(idx(k)).List = ones(1,prod(Tree.Val(arr)));
%     hw2PrintFactor(Mu(idx(k)),NameList);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calibrate the tree and compute the marginals for each variable
[Beta, Mu] = hw2CalibrateCliqueTree(Tree, Factors);
disp( 'Printout all beliefs Beta_i' );
for i=1:length(Beta)
    hw2PrintFactor(Beta(i),NameList);
end

% %Calculate the marginals
% %Manually
% P_1 = hw2FactorMarginalize( Beta(1), [2] );
% P_2 = hw2FactorMarginalize( Beta(1), [1] );
% P_3 = hw2FactorMarginalize( Beta(2), [2 4] );
% P_4 = hw2FactorMarginalize( Beta(2), [2 3] );
% P_5 = hw2FactorMarginalize( Beta(4), [3] );
% P_6 = hw2FactorMarginalize( Beta(3), [3 4] );
% P_7 = hw2FactorMarginalize( Beta(5), [4 6] );
% P_8 = hw2FactorMarginalize( Beta(6), [6 7] );
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
    hw2PrintFactorNormalized(P(i),NameList);
end
