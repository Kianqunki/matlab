function P = simpson(f,a,b,n)
%f=name of function, a=start value, b=end value, n=number of
%iterations

% Setting up
h=(b-a)/n;

S=f(a);
i=1:2:n-1;
x=a+h.*i;
y=f(x);

% Start computing
S=S+4*sum(y);
i=2:2:n-2;
x=a+h.*i;
y=f(x);
S=S+2*sum(y);
S=S+f(b);

P=h*S/3;

end
