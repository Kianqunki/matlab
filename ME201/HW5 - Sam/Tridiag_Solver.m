function x=Tridiag_Solver(e,f,g,b);

    N=length(f);
    
    for k=2:N
        
        mult=e(k)/f(k-1);
        
        f(k)=f(k)-mult*g(k-1);
        
        b(k)=b(k)-mult*b(k-1);
        
    end
    
    x(N)=b(N)/f(N);
    
    for k=N-1:-1:1;
        
        x(k)=(b(k)-g(k)*x(k+1))/f(k);
        
    end
        
  % End of program