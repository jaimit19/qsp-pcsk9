
clear all

% Load baseline model & get configuration
sbioloadproject('antiPCSK9_gadkar.sbproj', 'm1') ;
cs = getconfigset(m1);

% selection of dose group
DoseGroup_Index=5;  % interger number from 1-12 (see list below)

% Index:        Name:                 
%    1         10mg anti-PCSK9        
%    2         40mg anti-PCSK9       
%    3         150mg anti-PCSK9       
%    4         300mg anti-PCSK9     
%    5         400mg anti-PCSK9     
%    6         600mg anti-PCSK9       
%    7         800mg anti-PCSK9       
%    8         400mgQ4W anti-PCSK9   
%    9         400mgQ8W anti-PCSK9    
%    10        800mgQ8W anti-PCSK9     
%    11        200mgQ8W anti-PCSK9   
%    12        800mgQ12W anti-PCSK9   
               
% simulation of model
SimTime=100;  % This is the simulation end time
cs.StopTime=SimTime;
set(cs.SolverOptions, 'OutputTimes',0:SimTime)
DoseVar=m1.Dose(DoseGroup_Index);
simData = sbiosimulate(m1, cs,DoseVar);
[T1,X1] = selectbyname(simData, {'total_antipcsk9','LDLp'});
 
% plot simulation results
header=strcat(m1.dose(DoseGroup_Index).Name,' dose');
figure();
subplot(2,1,1);
semilogy(T1,X1(:,1),'LineWidth',2);
xlabel('Time (days)','FontSize',14)
ylabel('Total aPCSK9 (\mug/mL)','FontSize',14)
set(gca,'FontSize',12)
title(header,'FontSize',16);
subplot(2,1,2);
plot(T1,X1(:,2),'LineWidth',2);
xlabel('Time (days)','FontSize',14)
ylabel('LDLc (% of baseline)','FontSize',14)
set(gca,'FontSize',12)
