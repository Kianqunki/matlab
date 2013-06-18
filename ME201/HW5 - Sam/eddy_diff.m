function K=eddy_diff(ustar,L,z)

% 1st Approach
if L<=0;
    
    K=1.35*0.35*ustar*z/((1-15*z/L)^-0.25);
else
    K=0.35*ustar*z/(0.74+4.7*z/L);
end

% END OF PROGRAM
