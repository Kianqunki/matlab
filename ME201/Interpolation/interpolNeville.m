function yout = interpolNeville(x, y, xin)
% Input: 
% x, y: table of values. N-by-1 vectors.
% xin: the value at which the function is interpolated
% Output:
% yout: the interpolated value

% Optional error checking: 
% 1) x, y has the same length? More than some minimum number?
% 2) No identical values.

% find the index of value in the table closest to xin
[diff, idx] = min(abs(x-xin));
N = length(x);

% initialize
C = y;
D = y;
yout = y(idx);

% starting to compute C and D iteratively
for m=1:N-1
    for i=1:N-m
        if (x(i)-x(i+m) == 0)
            disp( 'ERROR: Identical values in the value tables' );
            return;
        else
            factor = (C(i+1) - D(i))/(x(i)-x(i+m));
            C(i) = factor*(x(i)-xin);
            D(i) = factor*(x(i+m)-xin);
        end
    end
    
    % decide which error to add into yout
    % Rule: branch such that it in the middle of interpolation window
    if ( idx < (N-m)/2 || idx == 1 )
        yout = yout + C(idx);
    else
        idx = idx - 1;
        yout = yout + D(idx);        
    end
end

end