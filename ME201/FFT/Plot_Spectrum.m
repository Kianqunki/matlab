function Plot_Spectrum(y,Fs)

    n=length(y);

    Y=y/n;
    
    AbsY=2*abs(Y(1:n/2+1));
    
    freq=Fs/2*linspace(0,1,n/2+1);
    
    % Plot single sided spectrum
    
    h=figure;
    
    plot(freq,AbsY);
    
    title('Amplitude spectrum of y(t)');
    
    xlabel('Frequency (Hz)');
    
    ylabel('|Y(f)|');
    
    grid on;
    

    
    