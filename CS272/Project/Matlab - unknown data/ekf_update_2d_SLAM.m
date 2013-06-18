function [X_k1k1, P_k1k1, descriptors, dict] = ekf_update_associate_2d_SLAM(X_k1k, P_k1k, z, R, descriptors, dict)

% divide into two steps
% determine which landmarks are new, which ones are old
% update with the old landmarks first
% Compute C^T(phi_R), equation 4.36
J = [0 -1; 1 0];
X_k1k1 = X_k1k;
P_k1k1 = P_k1k;
T1=chi2inv(0.95,2);
T2=chi2inv(0.99999,2);

% for each measurements
for i=1:size(z.measurements,2)
	%     if (isempty(descriptors)) % no landmark to associate
	%         % just initialize
	%         % TODO: all of first batch or just first reading?
	%         C = [cos(X_k1k1(3)) -sin(X_k1k1(3));
	%             sin(X_k1k1(3)) cos(X_k1k1(3))];
	%         z_k1 = z.measurements(:,i);
	%         X_L = X_k1k1(1:2)+C*z_k1;
	%         HR_k1 = [ -C' -C'*J*(X_L - X_k1k1(1:2))];
	%         HL_k1 = C';
	%
	%         entryNum = size(X_k1k1,1);
	%         X_k1k1 = [X_k1k1; X_L];  %equation 4.69, 4.70
	%         Prr_k1 = P_k1k1(1:3,1:3);
	%         P_k1k1_new = zeros(entryNum+2);
	%         P_k1k1_new(1:entryNum,1:entryNum) = P_k1k1;
	%         P_k1k1_new(entryNum+1:entryNum+2,1:entryNum) = -HL_k1'*HR_k1*P_k1k1(1:3,1:entryNum);
	%         P_k1k1_new(1:entryNum,entryNum+1:entryNum+2)=-P_k1k1(1:entryNum,1:3)*HR_k1'*HL_k1;
	%         P_k1k1_new(entryNum+1:entryNum+2,entryNum+1:entryNum+2) = HL_k1'*(HR_k1*Prr_k1*HR_k1'+R)*HL_k1;
	%
	%         %update the output
	%         P_k1k1 = P_k1k1_new;
	%         landmarkName = length(descriptors)+1;
	%         descriptors = [descriptors landmarkName];
	%
	%         %DEBUG
	%         dict=[dict z.descriptors(i)];
	%     else
	% for all discovered landmarks
	check=zeros(1,length(descriptors));
	for j=1:length(descriptors)
		lmIdx=descriptors(j)
        X_k1k1
        input('Wait')

		% put C here for clarity purpose
		C = [cos(X_k1k1(3)) -sin(X_k1k1(3));
			sin(X_k1k1(3)) cos(X_k1k1(3))];
		X_L = X_k1k1(2*lmIdx+2:2*lmIdx+3);   % estimated location of landmark
		% first three is rover location, every next two is landmark'
		% location.

		z_k1 = z.measurements(:,i);
		z_hat_k1 = C'*(X_L - X_k1k1(1:2));
		r = z_k1 - z_hat_k1;
		HR_k1 = [ -C' -C'*J*(X_L - X_k1k1(1:2))];
		HL_k1 = C';
		H_k1 = zeros(2,length(X_k1k1));
		H_k1(:,1:3) = HR_k1;
		H_k1(:,2*lmIdx+2:2*lmIdx+3) = HL_k1;

		S_k = H_k1*P_k1k1*H_k1' + R;    % this can be faster
		check(j)=r'*inv(S_k)*r;    %r'*inv(S_k)*r
	end

	idx = find(check>T2);
	if (length(idx)==length(descriptors))% | length(descriptors)==0)
		% initialize
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
		landmarkName = length(descriptors)+1;
		descriptors = [descriptors landmarkName];

		%DEBUG
		dict=[dict z.descriptors(i)];
	end

	idx = find(check<T1);
	if (length(idx)==1)
		% associate
		lmIdx = descriptors(idx);
		% DEBUG
		fprintf(1,'Landmark ID. Associate with:%d Truth:%d.\n', z.descriptors(i), dict(lmIdx));
		if (z.descriptors(i) ~= dict(lmIdx) )
			reply = input('Wrong guess. Press any key to continue', 's');
		end

		C = [cos(X_k1k1(3)) -sin(X_k1k1(3));
			sin(X_k1k1(3)) cos(X_k1k1(3))];
		X_L = X_k1k1(2*lmIdx+2:2*lmIdx+3);   % estimated location of landmark
		% first three is rover location, every next two is landmark'
		% location.

		z_k1 = z.measurements(:,i);
		z_hat_k1 = C'*(X_L - X_k1k1(1:2));
		r = z_k1 - z_hat_k1;
		HR_k1 = [ -C' -C'*J*(X_L - X_k1k1(1:2))];
		HL_k1 = C';
		H_k1 = zeros(2,length(X_k1k1));
		H_k1(:,1:3) = HR_k1;
		H_k1(:,2*lmIdx+2:2*lmIdx+3) = HL_k1;

		S_k = H_k1*P_k1k1*H_k1' + R;
		%K_k1 = P_k1k1*H_k1'*inv(S_k);
		K_k1 = P_k1k1*H_k1'/(S_k);
		X_k1k1 = X_k1k1 + K_k1*r;
		P_k1k1 = P_k1k1 - K_k1*H_k1*P_k1k1;
	end
	%else do nothing
	%     end
end

end