% % BIOM1010 Tutorial: Excitable Tisue Modelling (Mechanisms)
% % Quantitative Descriptions of Neuronal Membrane Potential 
% % Tianruo Guo & Nigel Lovell 7/9/2017

clear
clc
close all

global Stim_Amp Stim_Onset Stim_Offset 

% % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % %
Stim_Amp_all =[18] ;  % stimulus amplitude (pA/cm^2)
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

[time,Out] = ode15s('Mechanisms_function',[t0:mid:tend],Initial);  
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
 
GNa = gNa*(m.^3).*h/0.00002/100000000;
GK = gK*(n.^4)/0.00002/100000000;             % conductance for sodium and potassium channels
% % % Assume the opening conductance 
% % % for sodium and potassium channels are 20 pS
% % % 20 pS =0.00002 uS 
% % % 1 cm^2=100000000 um^2
% % % gNa gK with Unit of uS/cm^2   
%plot of membrane potential and stimulus
 
subplot(2,1,1);
plot(time,V);hold on;
ylim([-90 50]);
xlabel('time (s)')
ylabel('Membrane Potential (mV)')

subplot(2,1,2);
plot(time,GNa, time, GK);hold on;
% ylim([-90 50]);
xlabel('time (s)')
ylabel('Number of opening ion channels/ um^2 ')
legend('sodium channel', 'potassium channel')

end

