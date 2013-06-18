xo = [0 1]';
sigma_p = 1;
sigma = 1;
zi = 1;
H = [-1 1];

figure
N = 1;
C = sigma_p^2*( eye(2) - H'*H*1/(2+sigma^2/N/sigma_p^2) );
x = C*(1/sigma_p^2*xo + 1/sigma^2*N*zi*H');
plot_2D_ellipse( x, C, 0.95, 'r' );


N = 2
C = sigma_p^2*( eye(2) - H'*H*1/(2+sigma^2/N/sigma_p^2) );
x = C*(1/sigma_p^2*xo + 1/sigma^2*N*zi*H');
plot_2D_ellipse( x, C, 0.95, 'g' );


N = 5
C = sigma_p^2*( eye(2) - H'*H*1/(2+sigma^2/N/sigma_p^2) );
x = C*(1/sigma_p^2*xo + 1/sigma^2*N*zi*H');
plot_2D_ellipse( x, C, 0.95, 'b' );


N = 100
C = sigma_p^2*( eye(2) - H'*H*1/(2+sigma^2/N/sigma_p^2) );
x = C*(1/sigma_p^2*xo + 1/sigma^2*N*zi*H');
plot_2D_ellipse( x, C, 0.95, 'k' );