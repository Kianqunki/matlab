function particles= pf_resample( particles, Nmin )
% Resample particles in the particle filter

particleNum = length(particles);
% compute the normalized importance weights
weights = [particles.weight];
weights = weights/sum(weights);
% compute the effective number of particles
Neff = 1/sum(weights.*weights);

if (Neff < Nmin)
    % Stochastic universal resampling: Algorithm 2
    % construct CDF: cumweights = c in the paper
    cumweights = cumsum(weights);
    % generate N random numbers: 
    div = 0:1/particleNum:(1-1/particleNum);
    s = rand(1)/particleNum + div; % u in the paper
%     s = rand(1,particleNum)/particleNum + div;

    % index of particles to be kept/parents
    keep = zeros(1,particleNum);
    
    i=1;
    for j=1:particleNum
        while ( i <= particleNum & s(j) > cumweights(i) )
            i = i+1;
        end
        
        % assign parent
        keep(j) = i;
    end
    
    % assign samples
    particles = particles(keep);
    % assign weights
    for i=1:particleNum
        % reset the weights
        particles(i).weight = 1/particleNum;
    end
else
    % just normalize the weights
    for i=1:particleNum
        % reset the weights
        particles(i).weight = weights(i);
    end
end
