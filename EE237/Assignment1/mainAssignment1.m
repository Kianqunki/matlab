clear all
close all

% Construct X = [q1, q2, q1', q2']
% Then, we have Xo = [0, 0, 0, 0]
Xo = [0, 0, 0, 0]';


% Use the ODE solver
options = odeset( 'RelTol', 1e-6 );
[t_out, X_out] = ode45(@DynamicFunction, [0 100], Xo, options);


% Draw animation
% lengths of the links
l1 = 0.8;
l2 = 0.6;

figure(2)
axis on, axis equal
axis([-0.5 1.5 -0.5 1.5])
hold on

for i=1:length(X_out)
    % top of link 1 position wrt origin
    G_p_l1 = [l1*cos(X_out(i,1));
            l1*sin(X_out(i,1))];
    % top of link 2 wrt top of link 1
    l1_p_l2 = [l2*cos(X_out(i,1) + X_out(i,2));
                l2*sin(X_out(i,1) + X_out(i,2))];
    % top of link 2 wrt origin
    G_p_l2 = G_p_l1 + l1_p_l2;
    
    links = [[0 0]', G_p_l1, G_p_l2];
    
    % clear the plot
    cla;
    % new plot
    hold on
    
    line( links(1,:), links(2,:) )
    plot( links(1,:), links(2,:), 'ro' )
    
    hold off
    drawnow
    
end

figure(1)
hold on
plot( t_out, X_out(:,1) )
plot( t_out, X_out(:,2), 'r' )
legend( 'q1', 'q2' )
hold off

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%