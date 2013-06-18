function Diffusion(h,k,alpha,L,Tm_limit,frtim,Ta,To)

    delx=0.1*L; Fo=0.1;
    
    delt=Fo*delx^2/alpha;
    
    Bi=h*delx/k;
    
    N=round(L/delx)+1;
    
    Nt=round(Tm_limit/delt)+1;
    
    Nprint=round(frtim*Tm_limit/delt)+1;
    
    x(1)=0.0; tm(1)=0.0;
    
    for ix=1:N-1;
      
       x(ix+1)=x(ix)+delx;
       
    end 
    
    for it=1:Nt-1;
      
       tm(it+1)=tm(it)+delt;
       
    end  
    
    Nb=N-2;
    
    e(1:Nb)=-Fo; g(1:Nb)=-Fo;
    
    f(1:Nb)=(1+2*Fo); 
    
    f(1)=f(1)-Fo/(1+Bi);
    
    Temp(1:N,1)=To-Ta;
    
    for it=1:Nt-1
        
        b=Temp(2:N-1,it);
        
        b(Nb)=b(Nb)+Fo*(To-Ta);
        
        Temp(2:N-1,it+1)=Tridiag_Solver(e,f,g,b);
        
        Temp(1,it+1)=Temp(2,it+1)/(1+Bi);
        
        Temp(N,it+1)=To-Ta;
        
    end
    
    Temp=Temp+Ta;
    
    % Analytical solution
    
    for it=1:Nt
        
        for ix=1:N
            
            aa=h/k; bb=sqrt(alpha*tm(it));
            
            cc=x(ix)/(2*bb+1.0e-08);
            
            Tanal(ix,it)=erf(cc)+...
                exp(aa*x(ix)+(aa*bb)^2)*(1-erf(cc+aa*bb));
            
        end
        
    end
    
    Tanal=(To-Ta)*Tanal+Ta;
    
    
    % Plot results
    
    h=figure;
    
    plot(x,Temp(:,Nprint),'-r');
    
    hold on;
    
    plot(x,Tanal(:,Nprint),'-b');
    
    grid on;
    
    xlabel('Distance'); ylabel('Temperature');
    
    legend('Numerical', 'Analytical',2);
    
 % End of program
    
    
    
    
    
    