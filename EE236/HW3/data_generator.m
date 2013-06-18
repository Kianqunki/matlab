function [U, X] = data_generator( params, N )

% dimension of observation and input
dim = 1;

% generate N bipolar, white, pseudorandom, unit magnitude.
U = rand(dim,N);
for i=1:N
    if (U(dim,i) > 0.5 )
        U(dim,i) = 1;
    else
        U(dim,i) = -1;
    end
end

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