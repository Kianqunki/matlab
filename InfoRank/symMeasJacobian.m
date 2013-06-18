function out = symMeasJacobian( i, j, S, r, l, params )

lmDim = params.lmDim;
rpDim = params.rpDim;
J = params.J;
K = params.K;
N = params.N;

out = sym(zeros(lmDim, rpDim*K + lmDim*N));

if S(i,j) == 0
    return;
else
    temp1 = -sym(eye(2));
    temp2 = -J*( r(1:2,i) - l(1:2,j) );
    temp = [temp1 temp2];
    H_L = [cos(r(3,i)) sin(r(3,i)); -sin(r(3,i)) cos(r(3,i))]; % H_L = C'(phi) only
    H_R = H_L*temp;
    out(:,(i-1)*rpDim+1:i*rpDim) = H_R;
    out(:,rpDim*K + (j-1)*lmDim+1:rpDim*K+j*lmDim) = H_L;
end

end