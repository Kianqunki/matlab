function y=Generate_Signal(Fs,Length,a1,f1,a2,f2,noise,Nplot);

    delt=1/Fs;
    
    t=(0:Length-1)*delt;
    
    % Sum of two frequencies
    
    x=a1*sin(2*pi*f1*t)+a2*sin(2*pi*f2*t);
    
    y=x+noise*randn(size(t));
    
    % Plot signal;
    
    h=figure;
    
    plot(t(1:Nplot),y(1:Nplot));
    
    title('Signal corrupted with zero-mean random noise');
    
    xlabel('time in secs');
    
    grid on;

    