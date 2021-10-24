
function Y = Mechanisms_function(t,Yin)
global Stim_Amp Stim_Onset Stim_Offset 

Y = zeros(4,1);                 %initialize output
 
V = Yin(1);                     %assigning variable to input
n = Yin(2);
m = Yin(3);
h = Yin(4);
 
gNa = 120000 ;                   %parameter
gK = 36000;
gL = 300;
VNa = 55;
VK = -72;
VL = -49;
Cm = 1;
 


if (t >= Stim_Onset) && (t < Stim_Offset)  %stimulus
    Istim = -Stim_Amp*1000;
else
    Istim = 0;
end
 
an = 10*(V+50)/(1-exp(-(V+50)/10));         %gating variables
bn = 125*exp(-(V+60)/80);
am = 100*(V+35)/(1-exp(-(V+35)/10));
bm = 4000*exp(-(V+60)/18);
ah = 70*exp(-(V+60)/20);
bh = 1000/(1+exp(-(V+30)/10));
 
%ODE for action potential,n-variable, m-variable and L-variable
%respectively
Y(1) = -(1/Cm)*((gNa*(m^3)*h*(V-VNa))+(gK*(n^4)*(V-VK))+(gL*(V-VL))+Istim);
Y(2) = an*(1-n)-bn*n;
Y(3) = am*(1-m)-bm*m;
Y(4) = ah*(1-h)-bh*h;end
