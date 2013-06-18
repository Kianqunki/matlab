function out = symPropJacobian( i, S, r, params )

lmDim = params.lmDim;
rpDim = params.rpDim;
J = params.J;
K = params.K;
N = params.N;

out = sym(zeros(rpDim, rpDim*K + lmDim*N));

if i > size(S,1)-1
    return;
else
    % compute Phi_{R_k}
    temp2 = sym(zeros(3));
    temp2(1:2,3) = J*( r(1:2,i+1) - r(1:2,i) );
    temp2 = temp2 + eye(3);
    out(:,(i-1)*rpDim+1:i*rpDim) = -temp2;
    out(:,i*rpDim+1:(i+1)*rpDim) = eye(3);
end

end