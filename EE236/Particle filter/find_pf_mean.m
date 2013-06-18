function X_m = find_pf_mean( particles );

sum = zeros(size(particles(1).Xr));
for i = 1:length(particles)
	sum = sum + particles(i).weight*particles(i).Xr;
end

X_m = sum;
	
	