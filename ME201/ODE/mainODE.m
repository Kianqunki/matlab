clear all
close all

global R_CONSTANT

% part a
R_CONSTANT = 20;
options = odeset('RelTol',1e-3,'AbsTol',[1e-3 1e-3 1e-3]);
[T,X] = ode45(@hw3Lorenz,[0:0.1:25],[1 1 1],options);
% the trajectories in in xy, xz, and yz planes
figure(11)
plot(X(:,1),X(:,2),'b*')
xlabel ('x')
ylabel ('y')
title ( 'Question a, r = 20' );
figure(12)
plot(X(:,1),X(:,3),'r*')
xlabel ('x')
ylabel ('z')
title ( 'Question a, r = 20' );
figure(13)
plot(X(:,2),X(:,3),'g*')
xlabel ('y')
ylabel ('z')
title ( 'Question a, r = 20' );

% plot the outputs x, y, z versus t
figure(14)
plot(T,X(:,1),'b*')
hold on
plot(T,X(:,2),'r*')
plot(T,X(:,3),'g*')
xlabel ('t')
legend( 'x', 'y', 'z');
title ( 'Question a, r = 20' );

% part b
R_CONSTANT = 28;
options = odeset('RelTol',1e-3,'AbsTol',[1e-3 1e-3 1e-3]);
[T,X] = ode45(@hw3Lorenz,[0:0.1:25],[1 1 1],options);
% the trajectories in in xy, xz, and yz planes
figure(21)
plot(X(:,1),X(:,2),'b*')
xlabel ('x')
ylabel ('y')
title ( 'Question b, r = 28' );
figure(22)
plot(X(:,1),X(:,3),'r*')
xlabel ('x')
ylabel ('z')
title ( 'Question b, r = 28' );
figure(23)
plot(X(:,2),X(:,3),'g*')
xlabel ('y')
ylabel ('z')
title ( 'Question b, r = 28' );

% plot the outputs x, y, z versus t
figure(24)
plot(T,X(:,1),'b*')
hold on
plot(T,X(:,2),'r*')
plot(T,X(:,3),'g*')
xlabel ('t')
legend( 'x', 'y', 'z');
title ( 'Question b, r = 28' );

% plot the 3D trajectory
figure(27)
plot3(X(:,3),X(:,2),X(:,1))
xlabel ('z')
ylabel ('y')
zlabel ('x')
title ( 'Question b, r = 28' );

% part c
R_CONSTANT = 28;
options = odeset('RelTol',1e-3,'AbsTol',[1e-3 1e-3 1e-3]);
[T,X] = ode45(@hw3Lorenz,[0:0.1:25],[6 6 6],options);
% % the trajectories in in xy, xz, and yz planes
% figure(31)
% plot(X(:,1),X(:,2),'b*')
% xlabel ('x')
% ylabel ('y')
% title( '2D plot for initial point at (6,6,6)');
% figure(32)
% plot(X(:,1),X(:,3),'r*')
% xlabel ('x')
% ylabel ('z')
% title( '2D plot for initial point at (6,6,6)');
% figure(33)
% plot(X(:,2),X(:,3),'g*')
% xlabel ('y')
% ylabel ('z')
% title( '2D plot for initial point at (6,6,6)');
% 
% % plot the 3D trajectory
% figure(34)
% plot3(X(:,3),X(:,1),X(:,2))
% xlabel ('z')
% ylabel ('x')
% zlabel ('y')
% title( '3D plot for initial point at (6,6,6)');

[T2,X2] = ode45(@hw3Lorenz,[0:0.1:25],[6 6.01 6],options);
% % the trajectories in in xy, xz, and yz planes
% figure(35)
% plot(X(:,1),X(:,2),'b*')
% xlabel ('x')
% ylabel ('y')
% title( '2D plot for initial point at (6,6.01,6)');
% figure(36)
% plot(X(:,1),X(:,3),'r*')
% xlabel ('x')
% ylabel ('z')
% title( '2D plot for initial point at (6,6.01,6)');
% figure(37)
% plot(X(:,2),X(:,3),'g*')
% xlabel ('y')
% ylabel ('z')
% title( '2D plot for initial point at (6,6.01,6)');

% plot the 3D trajectory
figure(38)
plot3(X(:,3),X(:,2),X(:,1))
hold on
plot3(X2(:,3),X2(:,2),X2(:,1), 'r')
xlabel ('z')
ylabel ('y')
zlabel ('x')
legend( '(6,6,6)', '(6,6.01,6)');
title( '3D plot for initial points at (6,6,6) and (6,6.01,6)');