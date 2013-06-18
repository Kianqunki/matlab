function [deltaP, covP] = estimate_BLUE( Antenna, c, noise_std, t, P )

% number of data
num = size(Antenna,2);

% compute R_minus, tau, and Jacobian h given A and P
R_minus = zeros(num,1);
tau = zeros(num,1);
h = zeros(num,2);

for i=1:num
    rel_dist = Antenna(:,i) - P;
    R_minus(i) = norm(rel_dist);
    h(i,:) = -rel_dist'./R_minus(i);
    tau(i) = t(i) - R_minus(i)/c;
end

% compute z and H given tau(1:end)
z = zeros(num-1,1);
H = zeros(num-1,2);
for i=1:length(z)
    z(i) = tau(i+1) - tau(1);
    H(i,:) = (h(i+1,:) - h(1,:))./c;
end

% the matrix A in equation 6.27
A = [-1 1 0,
    -1 0 1];

% compute theta = deltaP
C_temp = inv(A*A');
deltaP = inv(H'*C_temp*H)*H'*C_temp*z;
covP = noise_std^2*inv(H'*C_temp*H);

end