function y=fftcr(x);

    x=x(:);
    
    n=length(x);
    
    omega=exp(-2*pi*i/n);
    
    j=0:n-1;
        
    k=j';
        
    F=omega.^(k*j);
        
    y=F*x;
        
   
    
    
    
    
        
        