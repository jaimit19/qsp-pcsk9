tic;
toc;

% Parameters of the antiPCKS9 model
function p =  antiPCSK9modelPars()
p.LDLparticleCE  =  0.92;      % milligram/nanomole 
p.circ_volume  = 5;         % liter              
p.clearance_hepatic_fraction  = 0.8;                          
p.baseline_hepatic_cholesterol   = 6000; %      milligram   ??       
p.maxSREBP2level  = 2 ;                            % ??
p.minSREBP2level = 0;                             % ??
p.LDLcClearanceRate =0.187; %1/day              
p.deciliter_to_liter = 10;       %deciliter/liter    
p.dilipidemic_index = 1;                            
p.PK_Ka  =  0.247; %     1/day              
p.PK_Kel_F  = 0.0443; %    1/day              
p. PK_V2_F  =  5640; %      milliliter         
p.PK_ComplexClearanceRate = 0.182; %     1/day              
p.PK_kon  = 0.84 ;     %1/nanomole/day     
p.PK_koff  = 1.11; %      1/day              
p.PK_V3  = 3990; %      milliliter         
p.PK_Q   = 615; %       milliliter/day     
p.Baselinepcsk9  = 281.94; %    nanogram/milliliter
p.LDLrIndClearanceRate = 0.02;       

p.AbsorptionFraction = 0.5;
p.CholesterolIntakeDiet = 300;
p.StatinEffectOnCholesterolSynthesis = 1;
p.BaselineCholesterolSynthesisRate = 800;
p.LDLparticleProdRate = 0.03;
p.LossRate = 0.23;
p.HDLcClearanceRate = 0.3;
p.pcsk9SynthesisRate = 9.4488;
p.pcsk9_synthesis.Vm_up = 2;
p.pcsk9_synthesis.Km_up = 1.5;
p.pcsk9_synthesis.Vm_down = 0.7;
p.pcsk9_synthesis.Km_down = 0.5;
p.pcsk9ClearanceRate = 2.48;
p.LDLr0 = 1;
p.gamma = 0;
p.LDLrSynthesis = 0.5;
p.LDLr_expression.Vm_up = 3;
p.LDLr_expression.Km_up = 1.5;
p.LDLr_expression.Vm_down = 0.7;
p.LDLr_expression.Km_down = 0.5;
p.LDLrClearance = 0.5;
p.pcsk9_on_LDLr = 0.75;
p.pcsk9_on_LDLr_range = 650;
p.antiPCSK9 = 1;
p.SREBP2 = 1;
p.circ_pcsk9_ngperml = 281.94;
%      LDLch = 80
%      antipcsk9_ugperml = 1.5e-05
%      complex_ugperml = 0
%      total_antipcsk9 = 1.5e-05
%      total_pcsk9 = 281.94
%      LDLp = 100
%      logPK = -4.8239
%      logTotalpcsk9 = 2.4502
end

function IC = getIntitialConditions()
hepatic_cholesterol = 6000;
LDLc = 4000;
circ_pcsk9 = 3.81;
surface_LDLr = 1;
antipcsk9 = 0.0001;
antipcsk9_dose = 0;
complex = 0;
peripheral = 0;
% HDLch = 50;
% LDLbaseline = 80;
IC = [hepatic_cholesterol, LDLc, circ_pcsk9, surface_LDLr, ...
    antipcsk9, antipcsk9_dose, complex, peripheral]; %...
    %HDLch, LDLbaseline];
end


function dydt = odefun(t,y,p)
dydt = zeros(length(y), 1);
hepatic_cholesterol = y(1);
LDLc = y(2);
circ_pcsk9 = y(3);
surface_LDLr = y (4);
antipcsk9 = y(5);
antipcsk9_dose = y(6);
complex = y(7) ;
peripheral = y(8);

% repeated assignment

SREBP2 = transform(hepatic_cholesterol,...
    cat(2, p.maxSREBP2level, p.minSREBP2level), ...
    p.baseline_hepatic_cholesterol, 3)';
    
% Fluxes
DietCholesterolAbsorption = p.AbsorptionFraction * p.CholesterolIntakeDiet;

CholesterolSynthesis = p.BaselineCholesterolSynthesisRate...
    * p.StatinEffectOnCholesterolSynthesis;

LDLFormation = p.LDLparticleProdRate * hepatic_cholesterol...
    * p.circ_volume * p.LDLparticleCE;

LDLcleranceToHepatic = p.LDLcClearanceRate * p.dilipidemic_index *...
    LDLc * surface_LDLr * p.clearance_hepatic_fraction;

CholesterolLost = p.LossRate * hepatic_cholesterol;

HDLclerance = p.HDLcClearanceRate * HDLch * p.circ_volume * p.deciliter_to_liter...
    * p.clearance_hepatic_fraction;

LDLcleranceToPeriphery = p.LDLcClearanceRate * p.dilipidemic_index * LDLc * (1- ...
    p.clearance_hepatic_fraction) *surface_LDLr;


pcsk9_synthesis = p.pcsk9SynthesisRate* SREBP2_reg(SREBP2,...
    p.pcsk9_synthesis.Vm_up,...
    p.pcsk9_synthesis.Km_up,...
    p.pcsk9_synthesis.Vm_down,...
    p.pcsk9_synthesis.Km_down);

pcsk9_clearance = p.pcsk9ClearanceRate * circ_pcsk9 ...
    * (surface_LDLr / p.LDLr0)^p.gamma;

%??
LDLr_expression = p.LDLrSynthesis * SREBP2_reg(SREBP2, ...
    p.LDLr_expression.Vm_up, p.LDLr_expression.Km_up, ...
    p.LDLr_expression.Vm_down, p.LDLr_expression.Km_down);

LDLr_clearance = LDLrClearance * surface_LDLr * ...
    transform(p.circ_pcsk9_ngperml, cat(2,(1 - p.pcsk9_on_LDLr), (1 + p.pcsk9_on_LDLr)), ...
    p.Baselinepcsk9, p.pcsk9_on_LDLr_range,'lin');

absorption = ((p.PK_Ka * antipcsk9_dose / p.PK_V2_F) * (1000 / 150) ) ...
    * antiPCSK9;

binding = p.PK_kon * antipcsk9 * circ_pcsk9 - p.PK_koff * complex;

antipcsk9_clearance = p.PK_Kel_F * antipcsk9;

complex_clearance = p.PK_ComplexClearanceRate * complex;

dose_compartment_clearance = p.PK_Ka * antipcsk9_dose;

distribution = p.PK_Q / p.PK_V2_F * (antipcsk9 - peripheral);

redistribution = p.PK_Q / p.PK_V3 * (peripheral - antipcsk9);

LDLrIndependentCleranceToHepatic = p.LDLrIndClearanceRate * ...
    p.clearance_hepatic_fraction * LDLc;

LDLrIndependentCleranceToPeriphery = p.LDLrIndClearanceRate * ...
    (1 - p.clearance_hepatic_fraction) * LDLc;


% ODEs
dhepatic_cholesterol_dt = DietCholesterolAbsorption + ...
    CholesterolSynthesis - LDLFormation + LDLcleranceToHepatic - ...
    CholesterolLost + HDLclerance + LDLrIndependentCleranceToHepatic;

dLDLc_dt = LDLFormation - LDLcleranceToHepatic - ...
    LDLcleranceToPeriphery - LDLrIndependentCleranceToHepatic - ...
    LDLrIndependentCleranceToPeriphery;

dcirc_pcsk9_dt = pcsk9_synthesis - pcsk9_clearance - binding;

dsurface_LDLr_dt = 1/p.antiPCSK9 * (LDLr_expression - LDLr_clearance);

dantipcsk9_dt = absorption - binding - antipcsk9_clearance - distribution;

dantipcsk9_dose_dt = -dose_compartment_clearance;

dcomplex_dt = binding - complex_clearance;

dperipheral_dt = -redistribution;


dydt(1) = dhepatic_cholesterol_dt;
dydt(2) = dLDLc_dt;
dydt(3) = dcirc_pcsk9_dt;
dydt(4) = dsurface_LDLr_dt;
dydt(5) = dantipcsk9_dt;
dydt(6) = dantipcsk9_dose_dt;
dydt(7) = dcomplex_dt ;
dydt(8) = dperipheral_dt;

end


function [Yout] = SREBP2_reg(SREBP2,Vmax_up, Km_up, Vmax_down, Km_down)

if ( SREBP2 >= 1 )
    In = (SREBP2-1);
    Yout = 1 + (Vmax_up-1) * In / (In + (Km_up-1)); 
end

if ( SREBP2 < 1 )
    In = (1-SREBP2);
    Yout = 1 - Vmax_down * In / (In + Km_down); 
end
 
end


function [Yout] = transform(Xin,Yrange,Xic50,varargin)

nVarargs = length(varargin);
Xscale='log'; % default is log
Xrange=3;     % log-order if log; actual range if linear
if (nVarargs==1) 
    Xrange=varargin{1}; 
end
if (nVarargs==2)
    Xrange=varargin{1}; Xscale=varargin{2};
end

if strcmp(Xscale,'log')
    expT=3.97865*Xrange^-0.995; % This is fitted to give 99% output at the extreme of the range specified
    Yout=( Yrange(2) - Yrange(1) ) * (Xin^expT/(Xin^expT+Xic50^expT)) + Yrange(1);
end

 if strcmp(Xscale,'lin')
    Xmin=Xic50-Xrange/2;Xmax=Xic50+Xrange/2;
    expT=3.97865*2^-0.995; % optimized for range of -1 to 1;
    Xmod=10^( 2*(Xin-Xmin)/(Xmax-Xmin) -1 )  ; 
    Xmodic50=10^( 2*(Xic50-Xmin)/(Xmax-Xmin) -1 );  
    Yout=( Yrange(2) - Yrange(1) ) * (Xmod^expT/(Xmod^expT+Xmodic50^expT)) + Yrange(1);
 end
 
end