clear all
close all

m = [1 2]';
v = [sqrt(2)/2 sqrt(2)/2]';

P1 = eye(2);
P2 = P1 - 0.5*(v*v');
P3 = P1 - 0.8*(v*v');

figure
plot_2D_ellipse(m, P1, 0.9, 'r');
plot_2D_ellipse(m, P2, 0.9, 'g');
plot_2D_ellipse(m, P3, 0.9, 'b');
axis equal

% Observations
% When the two variables are independent, the p-percentile error is bounded
% by a circle.
% The more dependence between the variables, the more elongated the
% p-percentile error ellipse is.