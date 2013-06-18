function particles= fastslam_resample( particles, Nmin )
% Resample particles in the particle filter

particleNum = length(particles);
% compute the normalized importance weights
weights = [particles.weight];
weights = weights/sum(weights);
% compute the effective number of particles
Neff = 1/sum(weights.*weights);

if Neff < Nmin
    % Stochastic universal resampling
    % do resampling
    cumweights = cumsum(weights);
    % generate P random numbers
    div = 0:1/particleNum:(1-1/particleNum);
    %s = rand(1)/particleNum + div;
    s = rand(1,particleNum)/particleNum + div;
    % index of particles to be kept
    keep = zeros(1,particleNum);
    
    j=1;
    for i=1:particleNum
        while ( j <= particleNum & s(j) <= cumweights(i) )
            keep(j) = i;
            j = j+1;
        end
    end
    
    particles = particles(keep);
    for i=1:particleNum
        % reset the weights
        particles(i).weight = 1/particleNum;
    end
end
