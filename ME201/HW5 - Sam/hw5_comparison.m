% Comparison of computed and measured crosswind integrated concentration
% divided by the source strength fro three downwind distances 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Monin-Obukhov length%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
close all
% Loading Datas
C_observed=load('C_observed.txt');

L=load('Monin_obukhov.txt');

ustar=load('ustar.txt');

% Computing the Concentrations from L,ustar given 
for i=1:60
    
    C_computed(i,1:3)=surface_release(1,0.46,@velocity,@eddy_diff,800,1500,0.008,L(i),ustar(i));
    
end

% Plot

loglog(1000*C_computed(:,1),C_observed(:,1),'o','MarkerFaceColor','g')
hold on

loglog(1000*C_computed(:,2),C_observed(:,2),'o','MarkerFaceColor','b')
loglog(1000*C_computed(:,3),C_observed(:,3),'o','MarkerFaceColor','r')

hold on 

%lines of 50% and 200%
loglog(1:1e4,1:1e4);loglog(1:1e4,2*(1:1e4));loglog(1:1e4,0.5*(1:1e4))

ylim([1 1e3])

xlim([1 1e3])

grid on

Title('Comparison of Computed and Measured Crosswind integrated concentration')

xlabel('Computed Concentration (\chi/Q)*1e3')

ylabel('Observed Concentration (\chi/Q)*1e3')
