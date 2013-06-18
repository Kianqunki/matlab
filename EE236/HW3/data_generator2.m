function [U, X] = data_generator2( params, N )

% dimension of observation and input
dim = 1;

% generate sum of sinusoids.
A1 = 1;
w1 = 0.1;
A2 = 2;
w2 = 0.2;
A3 = 3;
w3 = 0.3;
t = 1:N;

U = A1.*sin(w1*t) + A2.*sin(w2*t) + A3.*sin(w3*t);

% generate process noise w
w = randn(dim,N);

% generate observations x
trueTheta = params.trueTheta;
X = zeros(dim,N);
% this is the easieast method. For large p, augment U with zeroes.
for i = 1:N
    X(dim,i) = X(dim,i) + U(dim,i)*trueTheta(1);
    
    if ( i > 1 )
        X(dim,i) = X(dim,i) + U(dim,i-1)*trueTheta(2);
    end
    
    if ( i > 2 )
        X(dim,i) = X(dim,i) + U(dim,i-2)*trueTheta(3);
    end
    
    if ( i > 3 )
        X(dim,i) = X(dim,i) + U(dim,i-3)*trueTheta(4);
    end
    
    if ( i > 4 )
        X(dim,i) = X(dim,i) + U(dim,i-4)*trueTheta(5);
    end
    
    % add the noise
    X(dim,i) = X(dim,i) + w(dim,i);
end

end