function [Beta, Mu] = hw2MessageCliqueTree(Beta, Mu, i, j, idx);
%hw2MessageCliqueTree Prepare BU message passing between two cliques

% %DEBUG
% fprintf(1,'%d %d %d\n', i, j, idx );
% NameList = {'c' 'd' 't' 'g' 'i' 's' 'l' 'j'};
% disp( 'BU Message' );

var = setdiff(Beta(i).Variable,Mu(idx).Variable);
sigma_i = hw2FactorMarginalize(Beta(i), var);

% hw2PrintFactor(sigma_i, NameList);
% hw2PrintFactor(Mu(idx), NameList);
delta_i = hw2FactorDivide(sigma_i, Mu(idx) );
% hw2PrintFactor(delta_i, NameList);

Beta(j) = hw2FactorProduct( Beta(j), delta_i);
Mu(idx) = sigma_i;

% %DEBUG
% if (j == 1)
%     input( 'Wait' );
% end

end