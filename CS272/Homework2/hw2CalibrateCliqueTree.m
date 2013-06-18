function [Beta, Mu] = hw2CalibrateCliqueTree(Tree, Factors);
%hw2CalibrateCliqueTree Belief update for clique tree calibration.

%DEBUG
NameList = {'c' 'd' 't' 'g' 'i' 's' 'l' 'j'};

% %Initialize Tree
% number of cliques
N = length(Tree.nodes);

% compute assignment of factors to cliques
alpha = zeros(length(Factors), 1);

for i = 1:length(Factors)
    for j = 1:N
        if ( isempty(setdiff(Factors(i).Variable,Tree.nodes{j})))
            % no difference between factor and the clique, aka factor in
            % clique
            % store the index of the clique to alpha
            alpha(i) = j;
            break; 
        end
    end
end

% initialize cliques potentials and messages
% Procedure Initialize-CTree in pg. 367
Beta = repmat(struct('List', [], 'Val', [], 'Variable', []), N, 1);
Mu = repmat(struct('List', [], 'Val', [], 'Variable', []), N, N);

for i = 1:N
%     factorIndex = find(alpha == i)';
%     if(length(factorIndex) <1), 
%         continue;
%     end
%     
%     Beta(i) = Factors( factorIndex(1) ); % assign to the first factor
%     for j = factorIndex(2:end)
%         Beta(i) = hw2FactorProduct( Beta(i),Factors(j) );
%     end

    Beta(i).Variable = Tree.nodes{i};
    Beta(i).Val = Tree.Val(Beta(i).Variable);
    Beta(i).List = ones(1,prod(Beta(i).Val));
    factorIndex = find(alpha == i)';
    if (length(factorIndex) <1)
        continue;
    end
    
    for j = factorIndex(1:end)
        Beta(i) = hw2FactorProduct( Beta(i),Factors(j) );
    end
end
% for k = 1:N
%     hw2PrintFactor(Beta(k),NameList);
% end

edgeIdx = find( Tree.edges == 1);
[iSub jSub] = ind2sub([N N], edgeIdx);
% set Mu to 1
for k=1:length(edgeIdx)
    [arr, idx1, idx2] = intersect(Tree.nodes{iSub(k)},Tree.nodes{jSub(k)});
    Mu(edgeIdx(k)).Variable = arr;
    Mu(edgeIdx(k)).Val = Tree.Val(arr);
    Mu(edgeIdx(k)).List = ones(1,prod(Tree.Val(arr)));
%     hw2PrintFactor(Mu(edgeIdx(k)),NameList);
end
% end of Initialize-CTree

flag=true;
SMALL = 0.00001;

% perform clique tree calibration
while (flag)
    % propagate the message before checking
    randomNum = ceil(2*length(edgeIdx)*rand);
    k = ceil(length(edgeIdx)*rand);
    %BU_message
    % choose a random messaging direction
%     hw2PrintFactor(Mu(edgeIdx(k)), NameList);
    if (randomNum > length(edgeIdx))
        % iSub sending, jSub receiving
        [Beta, Mu] = hw2MessageCliqueTree(Beta, Mu, iSub(k),jSub(k), edgeIdx(k));
    else
        % iSub receiving, jSub sending
        [Beta, Mu] = hw2MessageCliqueTree(Beta, Mu, jSub(k),iSub(k), edgeIdx(k));
    end
%     hw2PrintFactor(Mu(edgeIdx(k)), NameList);
%     y=input('Wait');
    
    % check if the tree is calibrated
    count = 0;
    for k=1:length(edgeIdx)
        var = setdiff(Beta(iSub(k)).Variable,Mu(edgeIdx(k)).Variable);
        sigma_i = hw2FactorMarginalize(Beta(iSub(k)), var);
        var = setdiff(Beta(jSub(k)).Variable,Mu(edgeIdx(k)).Variable);
        sigma_j = hw2FactorMarginalize(Beta(jSub(k)), var);
        
%         hw2PrintFactor(sigma_i, NameList);
%         hw2PrintFactor(Mu(edgeIdx(k)), NameList);
%         y=input('Wait');
        
        if ( isempty(find(sigma_i.List~=Mu(edgeIdx(k)).List)) && isempty(find(sigma_j.List~=Mu(edgeIdx(k)).List)))           
            count = count+1;
        end
    end
    
    if (count==length(edgeIdx))
        % the tree is calibrated
        flag=false;
    end

end

end