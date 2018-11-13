clc
clear
close all

time = 15;
IC = zeros(1,4);

global a1 a2 a3 a4 a5 a6 a7 a8;

IC(1) = 500; %0.2-3 um is the range for [Trimeric G Protein]
IC(2) = 1; %[By]
IC(3) = 1; %[Galpha/GTP]
IC(4) = 1; %[Galpha/GDP]

a1 = 0; %0.001-0.5 um for [GPCR]
a2 = 20; %0.01-0.3 nm for [GAPS]
a3 = 0.01; %kass 
a4 = 20; %kdiss (in sec^-1)
a5 = 500; %K2 
a6 = 30; %RC complex
a7 = 5; %khydro (in sec^-1)
a8 = 2; %K3 

[t, yd] = ode45('gpcrfx',[0, time],IC);
figure(1)
plot(t, yd)
ylabel('Concentration (nm)')
xlabel('time (sec)')
legend('TriGP','By','GaGTP','GaGDP')