function [X_k1k1, P_k1k1, descriptors] = ekf_update_2d_SLAM(X_k1k, P_k1k, z, R, descriptors)

% divide into two steps
% determine which landmarks are new, which ones are old
% update with the old landmarks first
% Compute C^T(phi_R), equation 4.36
J = [0 -1; 1 0];
X_k1k1 = X_k1k;
P_k1k1 = P_k1k;

% batch updating of the discovered landmarks
counter=0;  % old landmark counter
H_k1 = zeros(2*length(z.descriptors),length(X_k1k1));
r = zeros(2*length(z.descriptors),1);

for i = 1:length(z.descriptors)
    idx = find(descriptors==z.descriptors(i),1); %the second parameter 1 is for faster processing
    % if the landmark is already discovered
    if (length(idx)==1)
        counter=counter+1;
        C = [cos(X_k1k1(3)) -sin(X_k1k1(3)); 
            sin(X_k1k1(3)) cos(X_k1k1(3))];
        X_L = X_k1k1(2*idx+2:2*idx+3);   % estimated location of landmark
        % first three is rover location, every next two is landmark'
        % location.
        
        z_k1 = z.measurements(:,i);
        z_hat_k1 = C'*(X_L - X_k1k1(1:2));
        r(2*counter-1:2*counter,1) = z_k1 - z_hat_k1;
        
        HR_k1 = [ -C' -C'*J*(X_L - X_k1k1(1:2))];
        HL_k1 = C';
        H_k1(2*counter-1:2*counter,1:3) = HR_k1;
        H_k1(2*counter-1:2*counter,2*idx+2:2*idx+3) = HL_k1;        
    end
end

R_big = zeros(2*counter);
for i=1:2:2*counter
    R_big(i:i+1,i:i+1)=R;
end
H_k1 = H_k1(1:2*counter,:);
r = r(1:2*counter,:);

S_k = H_k1*P_k1k1*H_k1' + R_big;
%K_k1 = P_k1k1*H_k1'*inv(S_k);
K_k1 = P_k1k1*H_k1'/(S_k);
X_k1k1 = X_k1k1 + K_k1*r;
P_k1k1 = P_k1k1 - K_k1*H_k1*P_k1k1;

% sequential initialization
for i = 1:length(z.descriptors)
    idx = find(descriptors==z.descriptors(i),1); %the second parameter 1 is for faster processing
    if (isempty(idx))
        % if the landmark is not discovered
        %change dimensions of X_k1k1 and P_k1k1
        C = [cos(X_k1k1(3)) -sin(X_k1k1(3));
            sin(X_k1k1(3)) cos(X_k1k1(3))];
        z_k1 = z.measurements(:,i);
        X_L = X_k1k1(1:2)+C*z_k1;
        HR_k1 = [ -C' -C'*J*(X_L - X_k1k1(1:2))];
        HL_k1 = C';
        
        entryNum = size(X_k1k1,1);
        X_k1k1 = [X_k1k1; X_L];  %equation 4.69, 4.70
        Prr_k1 = P_k1k1(1:3,1:3);
        P_k1k1_new = zeros(entryNum+2);
        P_k1k1_new(1:entryNum,1:entryNum) = P_k1k1;
        P_k1k1_new(entryNum+1:entryNum+2,1:entryNum) = -HL_k1'*HR_k1*P_k1k1(1:3,1:entryNum);
        P_k1k1_new(1:entryNum,entryNum+1:entryNum+2)=-P_k1k1(1:entryNum,1:3)*HR_k1'*HL_k1;
        P_k1k1_new(entryNum+1:entryNum+2,entryNum+1:entryNum+2) = HL_k1'*(HR_k1*Prr_k1*HR_k1'+R)*HL_k1;
        
        %update the output
        P_k1k1 = P_k1k1_new;
        descriptors = [descriptors z.descriptors(i)];
    end
end