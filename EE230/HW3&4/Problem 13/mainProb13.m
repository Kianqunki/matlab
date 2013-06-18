clear all
close all

load LS_data.mat

figure
plot(y)
hold on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Linear model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 13
A = [ones(100,1) t'];
[Q, R] = qr(A,0);
% the parameters
a = R\(Q'*y')
linear_prediction = [1 120]*a

plot( A*a, 'r')

% Problem 14
linear_error_norm = norm(y'-A*a)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% quadratic model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 13
A = [ones(100,1) t' t'.^2];
[Q, R] = qr(A,0);
% the parameters
a = R\(Q'*y')
quadratic_prediction = [1 120 120^2]*a

plot( A*a, 'g')

% Problem 14
quadratic_error_norm = norm(y'-A*a)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cubic model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 13
A = [ones(100,1) t' t'.^2 t'.^3];
[Q, R] = qr(A,0);
% the parameters
a = R\(Q'*y')
cubic_prediction = [1 120 120^2 120^3]*a

plot( A*a, 'k')

% Problem 14
cubic_error_norm = norm(y'-A*a)