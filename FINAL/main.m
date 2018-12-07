%Main method that executes the olfactory system and then references another
%script to illustrate how an input of a specific odorant concentration
%would affect the elucidation of action potentials. This is done by
%referencing two scripts: olfModel.m and actionPotential.m 
clc; clear; close all;

%%Olfactory System
IC = zeros(1,12);
odorantsConc = [0.01 0.40 2 55 66.5 200 1000];
timeCell = cell(1, length(odorantsConc));
outputHolder = cell(1, length(odorantsConc));
drug = 1;

%Global variables are used as these terms remain constant across subsystems
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
a = 2.0*drug; %rate of AC activation by active Golf 
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
    [t, yd] = ode45('olfModel', [0:0.00001:time], IC);
    outputHolder{i} = yd;
    timeCell{i} = t;
end

%%Plotting output of second messenger concentrations from olfModel.m
figure(1)
for i = 1:length(outputHolder)
    plot(timeCell{i}, outputHolder{i}(:,1));
    hold on
end
xlabel('time (sec)')
ylabel('Active CNG channels (\muM)')
ylim([0 0.2])
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
ylim([0 4])
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
ylabel('Ca^{2+} concentration (\muM)')
ylim([0 0.35])
title('Concentration of active Ca^{2+} ions over time for varying odorant conc')
legend(['=' num2str(odorantsConc(1)) '\muM'], ['=' num2str(odorantsConc(2)) '\muM'], ...
    ['=' num2str(odorantsConc(3)) '\muM'], ['=' num2str(odorantsConc(4)) '\muM'], ...
    ['=' num2str(odorantsConc(5)) '\muM'], ['=' num2str(odorantsConc(6)) '\muM'], ...
    ['=' num2str(odorantsConc(7)) '\muM']);

%Plots the maximum values encountered from executing olfModel.m 
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
ylabel('Concentration (\muM)')
legend('CNG Channels', 'cAMP', 'Ca2+');
title('Maximum values of Conc vs Odorant Conc')
fprintf('go to AP \n')


figure(5)
for i = 1:4   
    subplot(2,2,i)
    Iapp = (outputHolder{i}(:,1) + outputHolder{i}(:,2))*80;
    [t, yd, gbar_k, gbar_na, m, h, n] = actionPotential(timeCell{end}(end)*1000, Iapp);
    freq(i) = freqCounter(yd, t, -50);
    plot(t, yd)
    xlabel('Time (ms)')
    ylabel('Membrane Voltage (mV)')
    ylim([-90 50])
end

figure(6)
for i = 5:7   
    subplot(2,2,(i-4))
    Iapp = (outputHolder{i}(:,1) + outputHolder{i}(:,2))*(80);
    [t, yd, gbar_k, gbar_na, m, h, n] = actionPotential(timeCell{end}(end)*1000, Iapp);
    freq(i) = freqCounter(yd, t, -50);
    plot(t, yd)
    xlabel('Time (ms)')
    ylabel('Membrane Voltage (mV)')
    ylim([-90 50])
end

%Plots conductances of ionic channels featured in actionPotential.m for
%specific odorant concentration. 
figure(7)
for i = 1:4
    subplot(2,2,i)
    Iapp = (outputHolder{i}(:,1) + outputHolder{i}(:,2))*(80);
    [t, yd, gbar_K, gbar_Na, m, h, n] = actionPotential(timeCell{end}(end)*1000, Iapp);
    p1 = plot(t,gbar_K*n.^4);
    hold on
    p2 = plot(t,gbar_Na*(m.^3).*h);
    legend([p1, p2], 'Conductance for Potassium', 'Conductance for Sodium')
    ylabel('Conductance')
    xlabel('time (ms)')
    title('Conductance for Potassium and Sodium Ions in Simulated Neuron')
end

figure(8)
for i = 5:7
    subplot(2,2,(i-4))
    Iapp = (outputHolder{i}(:,1) + outputHolder{i}(:,2))*(80);
    [t, yd, gbar_K, gbar_Na, m, h, n] = actionPotential(timeCell{end}(end)*1000, Iapp);
    p1 = plot(t,gbar_K*n.^4);
    hold on
    p2 = plot(t,gbar_Na*(m.^3).*h);
    legend([p1, p2], 'Conductance for Potassium', 'Conductance for Sodium')
    ylabel('Conductance')
    xlabel('time (ms)')
    title('Conductance for Potassium and Sodium Ions in Simulated Neuron')
end


%Plots sequences of saturations and changes in concentrations of second
%messengers against time (this assumes that we have a 10 second execution
%of the system as it plots the first 5 seconds, otherwise an error is
%thrown) 
figure(9)
ind = 1;
for i = [1 2 3 7]
    timeInd = timeCell{i} <= 5;
    times = timeCell{i}(timeInd);
    subplot(2,2,ind)
    for j = [6 1 9]
        plot(times,outputHolder{i}(timeInd,j));
        hold on
    end
    xlabel('Time (sec)')
    ylabel('Concentration (\muM)')
    title(['Concentrations vs Time when Oderant Concentration = ' num2str(odorantsConc(i)) '\muM'])
    if ind == 1
        legend({'cAMP','CNG Channels','Ca^{2+}'},'Location','northwest','AutoUpdate','off')
    else
        legend({'cAMP','CNG Channels','Ca^{2+}'},'AutoUpdate','off')
    end
    for j = [6 1 9]
        [pks,locs] = findpeaks(outputHolder{i}(:,j));
        plot([timeCell{i}(locs(1)),timeCell{i}(locs(1))],[0 pks(1)],'k:');
        hold on
    end
    ind = ind + 1;
end

%Subfunction that calculates the frequency values of action potentials for
%specific odorant concentrations. 
function output = freqCounter(Iapp, timeCell, threshold)
   pass = 0;
   startTime = 0;
   endTime = 0;
   width = 0;
   count = 1;
    for i = 1:length(Iapp)
       if Iapp(i) > threshold && pass == 0
           startTime = timeCell(i);
           pass = 1;
       end
       if Iapp(i) < threshold && pass == 1
           endTime = timeCell(i);
           pass = 0;
       end
       if startTime ~= 0 && endTime ~= 0
            width(count) = endTime - startTime;
            endTime = 0;
            startTime = 0;
            count = count + 1;
       end
    end
    width = 1000./width;
    output = mean(width);
end