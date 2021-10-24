% % BIOM1010 Tutorial: Excitable Tisue Modelling (Action Potentials) 
% % Quantitative Descriptions of Neuronal Membrane Potential 
% % Tianruo Guo & Nigel Lovell 7/9/2017

clear
clc
close all

global Stim_Amp Stim_Onset Stim_Offset

% % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % %
Stim_Amp_all =[1 2 3 4 6 8 12] ;  % stimulus amplitude (pA/cm^2)
Stim_Onset=0.002;  % Stimulus Onset (s)
Stim_Offset=0.004; % Stimulus offset (s)
% % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % %

Initial = [-60, 0.3177, 0.0529, 0.5961];           %initial condition and time
t0 = 0;
mid = 0.0001;
tend = 0.02;

for i=1:length(Stim_Amp_all)
    Stim_Amp=Stim_Amp_all(i);

[time,Out] = ode15s('ALL_or_NONE_function',[t0:mid:tend],Initial);  
V = Out(:,1);
n = Out(:,2);
m = Out(:,3);
h = Out(:,4);               %result for a neuron
Stimulus=zeros(length(V),1) ;
Stimulus(Stim_Onset/mid+1: Stim_Offset/mid+1)=Stim_Amp;

gNa = 120000 ;   % uS/cm^2
gK = 36000;      % uS/cm^2
gL = 0.3;        % uS/cm^2
VNa = 55;        % mV
VK = -72;        % mV
VL = -49;        % mV
Cm = 1;                     %parameters for currents
 
INa = (gNa*(m.^3).*h.*(V-VNa));
IK = (gK*(n.^4).*(V-VK));
IL = (gL*(V-VL));           %currents calculation for normal neuron
 
GNa = gNa*(m.^3).*h;
GK = gK*(n.^4);             % conductance for sodium and potassium channels

%plot of membrane potential and stimulus
 
subplot(2,1,1);
plot(time,Stimulus);hold on;
ylim([-2 18])
xlabel('time (s)')
ylabel('Stimulus(pA/cm^2)')
title('All or None Principle');

subplot(2,1,2);
plot(time,V);hold on;
ylim([-90 50]);
xlabel('time (s)')
ylabel('Membrane Potential (mV)')

% subplot(4,1,3);
% plot(time,GNa);hold on;
% % ylim([-90 50]);
% xlabel('time (s)')
% ylabel('Conductance ')
% 
% subplot(4,1,4);
% plot(time,GK);hold on;
% % ylim([-90 50]);
% xlabel('time (s)')
% ylabel('Conductance ')
end

