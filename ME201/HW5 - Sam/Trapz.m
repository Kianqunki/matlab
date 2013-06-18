%This Function calculates the numerical integration using Trapezoidal rule

function S=Trapz(f,h)

lenx=length(f);

S=h/2*(f(1)+f(lenx)+2*sum(f(2:lenx-1)));



    
    



