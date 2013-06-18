function plot_2D_ellipse(mean, covar, prob, color)

[U D V] = svd(covar);
mThres = chi2inv(prob,2);

% find the length of the major and minor axes
a = sqrt(mThres*D(1,1));
b = sqrt(mThres*D(2,2));

% generate all the angles for a full cycle
theta = [0:1/200:2*pi+1/200]; %1/20 is added to ensure full cycle generated

% generate an ellipse
state = zeros(2,length(theta));
state(1,:) = a*cos(theta);
state(2,:) = b*sin(theta);
% rotate it
X = V * state;
% translate it
X(1,:) = X(1,:) + mean(1);
X(2,:) = X(2,:) + mean(2);

plot(X(1,:),X(2,:),color);
hold on;