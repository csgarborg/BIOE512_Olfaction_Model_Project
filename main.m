clc; clear; close all;

%%Olfactory System
IC = zeros(1,12);
odorantsConc = [0.01 0.40 2  30 70 200 10000];
timeCell = cell(1, length(odorantsConc));
outputHolder = cell(1, length(odorantsConc));

global od rc os os2 pka AKmax ak2 PKmax pk2 go rg rg2 ac a ckk A3max A2max...
    ca apd PDmax cng CNmax cn2 cam cc cc2 CMmax ifa ef cac CLmax...
    cl2 cpd CPmax cp2 CDmax CKmax ck2 K1 K2 K3 N1 N2 N3 tempHold;

%Initializing global variables
rc = 1.0; %concentration of receptor
os = 0.1; %rate of odorant binding to receptor 
os2 = 2.0; %rate of odorant dissociation from receptor
pka = 1.0; %concentration of PKA 
AKmax = 0.1; %maximal rate of PKA activation by cAMP
ak2 = 1.0; %rate of PKA inactivation 
PKmax = 1.0; %maximal rate of receptor inactivation by active PKA
pk2 = 0.1; %rate of reversal of receptor inactivation 
go = 1.0; %concentration of Golf
rg = 10; %rate of activation of Golf by odorant—receptor complex
rg2 = 2.0; %rate of inactivation of Golf
ac = 5.0; %concentration of AC
a = 2.0; %rate of AC activation by active Golf 
ckk = 1.0; %concentration of CAM-kinase II
A2max = 10; %maximal rate of AC inactivation by CAM-kinase II 
A3max = 40; %maximal rate of AC inactivation by calcium ions
ca = 40; %rate of cAMP synthesis by active AC
apd = 1.0; %concentration of cAMP-PDE
PDmax = 10; %maximal rate of cAMP hydrolysis by cAMP-PDE
cng = 1.0; %concentration of CNG channels
CNmax = 5.0; %maximal rate of CNG channel activation by cAMP
cn2 = 10; %rate of CNG channel inactivation
cam = 1.0; %concentration of CAM
cc = 0.1; %rate of calcium binding to CAM
cc2 = 1.0; %rate of calcium dissociation from CAM
CMmax = 10; %maximal rate of CNG channel inactivation by Ca4-CAM
ifa = 20; %rate of calcium ion influx through CNG channels
ef = 10; %rate of calcium ion efflux by NCX 
cac = 1.0; %concentration of Ca—Cl channels
CLmax = 5.0; %maximal rate of Ca—Cl channel activation by calcium ions
cl2 = 10; %rate of Ca—Cl channel inactivation
cpd = 1.0; %concentration of CAM-PDE
CPmax = 10; %maximal rate of CAM-PDE activation by Ca4CAM 
cp2 = 0.01; %rate of CAM-PDE inactivation
CDmax = 20; %maximal rate of cAMP hydrolysis by CAM-PDE
CKmax = 10; %maximal rate of CAM-kinase II activation by Ca4CAM
ck2 = 0.01; %rate of CAM-kinase II inactivation
K1 = 3.4; %concentration of cAMP at which CNG channel activation is half-maximal 
K2 = 5.0; %concentration of calcium ions at which Ca—Cl channel activation is half-maximal 
K3 = 0.2; %concentration of calcium ions at which AC inhibition is half-maximal
N1 = 1.4; %Hill coefficients for the processes 
N2 = 2; %Hill coefficients for the processes 
N3 = 3; %Hill coefficients for the processes 

%Setting initial values
IC(1) = 0.0001; %x
IC(2) = 0.0001; %y
IC(3) = 0.001; %z
IC(4) = 0.001; %w
IC(5) = 0.001; %v
IC(6) = 0.001; %u
IC(7) = 0.001; %s
IC(8) = 0.001; %r
IC(9) = 0.001; %q
IC(10) = 0.001; %p
IC(11) = 0.001; %o
IC(12) = 0.001; %n

time = 10; %seconds
for i = 1:length(odorantsConc)
    od = odorantsConc(i); %concentration of odorants (VARIABLE ALONG WITH TIME EXPOSED)
    tempHold = od;
    [t, yd] = ode45('olfModel', [0:0.001:time], IC);
    outputHolder{i} = yd;
    timeCell{i} = t;
end

figure(1)
for i = 1:length(outputHolder)
    plot(timeCell{i}, outputHolder{i}(:,1));
    hold on
end
xlabel('time (sec)')
ylabel('Active CNG channels (\muM)')
title('Concentration of active CNG channels over time for varying odorant conc')
legend(['=' num2str(odorantsConc(1)) '\muM'], ['=' num2str(odorantsConc(2)) '\muM'], ...
    ['=' num2str(odorantsConc(3)) '\muM'], ['=' num2str(odorantsConc(4)) '\muM'], ...
    ['=' num2str(odorantsConc(5)) '\muM'], ['=' num2str(odorantsConc(6)) '\muM'], ...
    ['=' num2str(odorantsConc(7)) '\muM']);

figure(2)
for i = 1:length(outputHolder)
    plot(timeCell{i}, outputHolder{i}(:,6));
    hold on
end
xlabel('time (sec)')
ylabel('cAMP (\muM)')
title('Concentration of cAMP over time for varying odorant conc')
legend(['=' num2str(odorantsConc(1)) '\muM'], ['=' num2str(odorantsConc(2)) '\muM'], ...
    ['=' num2str(odorantsConc(3)) '\muM'], ['=' num2str(odorantsConc(4)) '\muM'], ...
    ['=' num2str(odorantsConc(5)) '\muM'], ['=' num2str(odorantsConc(6)) '\muM'], ...
    ['=' num2str(odorantsConc(7)) '\muM']);

figure(3)
for i = 1:length(outputHolder)
    plot(timeCell{i}, outputHolder{i}(:,9));
    hold on
end
xlabel('time (sec)')
ylabel('Ca2+ concentration (\muM)')
title('Concentration of active Ca2+ ions over time for varying odorant conc')
legend(['=' num2str(odorantsConc(1)) '\muM'], ['=' num2str(odorantsConc(2)) '\muM'], ...
    ['=' num2str(odorantsConc(3)) '\muM'], ['=' num2str(odorantsConc(4)) '\muM'], ...
    ['=' num2str(odorantsConc(5)) '\muM'], ['=' num2str(odorantsConc(6)) '\muM'], ...
    ['=' num2str(odorantsConc(7)) '\muM']);

figure(4)
maxValue = zeros(3, length(odorantsConc));
types = [1 6 9];
for i = 1:length(types)
    for j = 1:length(outputHolder)
        maxValue(i, j) = max(outputHolder{j}(:,types(i)));
    end
    semilogx(odorantsConc, maxValue(i, :), '*-');
    hold on
end
xlabel('Odorant Concentration (\mum)')
ylabel('Molecular concentration (\muM)')
title('Maximum values of Conc vs Odorant Conc')

%%Action Potential system

%Concentration of active ion channels multipled by a Ki that lies between a
%range of 8-300 uM from paper of following html: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1171109/pdf/000131.pdf
%on page 10/14. 
Iapp = (outputHolder{5}(:,1) + outputHolder{5}(:,2)*80);
    
[t, yd] = actionPotential(timeCell{end}(end)*100, timeCell{end}, Iapp);
figure(5)
plot(t, yd)
xlabel('Time (ms)')
ylabel('Membrane Voltage (mV)')
title('Action Potential graph over time of neurons ')