function [dX_dt] = DynamicFunction( t, X )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Parameters: from page 159
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% masses of the links
m1 = 1;
m2 = 0.8;
% lengths of the links
l1 = 0.8;
l2 = 0.6;
% distances of the centers of the links
l_c1 = 0.4;
l_c2 = 0.3;
% moments of inertia of the links
J1 = 0.0533;
J2 = 0.024;


% gravity
g = 9.8;
% Thetas
Thetas = [ m1*l_c1^2 + m2*(l1^2 + l_c2^2) + J1 + J2;
    m2*l1*l_c2;
    m2*l_c2^2 + J2;
    m1*l_c1 + m2*l1;
    m2*l_c2];

Kp = 5*eye(2);
Kd = 0.25*eye(2);

% And we have q_d = [30*pi/180, 60*pi/180]
Xd = [30*pi/180, 60*pi/180,0,0]';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parse the state vector for easy formulation
q1 = X(1);
q2 = X(2);
q1_d = X(3); % derivative of q1
q2_d = X(4); % derivative of q2


Mi = [Thetas(1) + 2*Thetas(2)*cos(q2),   Thetas(3) + Thetas(2)*cos(q2);
    Thetas(3) + Thetas(2)*cos(q2),     Thetas(3)                    ];
Ci = [-Thetas(2)*sin(q2)*q2_d       ,    -Thetas(2)*sin(q2)*(q1_d + q2_d)
    Thetas(2)*sin(q2)*q1_d       ,     0                         ];
gi = [Thetas(4)*g*cos(q1) + Thetas(5)*g*cos(q1+q2);
    Thetas(5)*g*cos(q1+q2)];

dX_dt = zeros(size(X));

% q_dot
dX_dt(1:2) = X(3:4);
% q_dot_dot
M_inv = inv(Mi);
dX_dt(3:4) = -M_inv*Ci*X(3:4) - M_inv*Kp*(X(1:2)-Xd(1:2)) - M_inv*Kd*X(3:4);

return