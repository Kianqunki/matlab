function x=Gauss_Elimination(A,b);

    [nx,ny]=size(A);
    
    aug=[A,b];  % Creates augmented matrix
    
    for k=1:nx-1
        
        [cm,im]=max(abs(aug(k:nx,k)));
        
        if cm~=0
        
            dum=aug(k,:);
        
            aug(k,:)=aug(im+k-1,:);
        
            aug(im+k-1,:)=dum;
        
            for i=k+1:nx
            
                mult=aug(i,k)/aug(k,k);
            
                for j=k:(nx+1)
               
                    aug(i,j)=aug(i,j)-mult*aug(k,j);
               
                end
        
            end
            
        end
        
    end
        
    AA=aug(:,1:nx); bb=aug(:,nx+1);
    
    x(nx)=bb(nx)/AA(nx,nx);
    
    for k=nx-1:-1:1
        
        sum=0.0;
        
        for j=(k+1):nx;
            
            sum=sum+AA(k,j)*x(j);
            
        end
        
        x(k)=(bb(k)-sum)/AA(k,k);
        
    end
    
  % End of program