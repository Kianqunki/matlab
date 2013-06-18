load QRmatrices.mat

[Q1, R1] = QR_Gram(A_good);
% Testing
Q1'*Q1 - eye(size(Q1,2))
A_good - Q1*R1

[Q2, R2] = QR_Householder(A_good);
% Testing
Q2'*Q2 - eye(size(Q2,2))
A_good - Q2*R2

% The two outputs are the same
Q1 - Q2
R1 - R2

pause

[Q1, R1] = QR_Gram(A_bad);
% Testing
Q1'*Q1 - eye(size(Q1,2))
A_bad - Q1*R1

[Q2, R2] = QR_Householder(A_bad);
% Testing
Q2'*Q2 - eye(size(Q2,2))
A_bad - Q2*R2