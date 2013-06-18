function F = trapezoid(f,a,b,n)
%f=name of function, a=start value, b=end value, n=number of 
%iterations

% Setting up
h=(b-a)./n;

S=f(a);
i=1:1:n-1;
x=a+h.*i;
y=f(x);

% Start computing
S=S+2.*sum(y);
S=S+f(b);

F=h.*S./2;

end