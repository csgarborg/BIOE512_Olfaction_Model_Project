clc; clear; close all;

%%Olfactory System
IC = zeros(1,12);
odorantsConc = [0.01 0.40 2 30 66.5 200 10000];
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

% %%Action Potential system
% 
% simulationTime = 1000; %in milliseconds
% deltaT=.01;
% t=0:deltaT:simulationTime;
% 
% 
% %===specify the external current I===
% changeTimes = 50; %in milliseconds
% currentLevels = 50; %Change this to see effect of different currents on voltage (Suggested values: 3, 20, 50, 1000)
% 
% %Set externally applied current across time
% %Here, first 500 timesteps are at current of 50, next 1500 timesteps at
% %current of zero (resets resting potential of neuron), and the rest of
% %timesteps are at constant current
% %  I(1:500) = currentLevels; I(501:2000) = 50; I(2001:numel(t)) = 50;
% Iapp = (outputHolder{i}(:,1) + outputHolder{i}(:,2))*80;
% % I(1:numel(t)) = 50;
% %Comment out the above line and uncomment the line below for constant current, and observe effects on voltage timecourse
% %I(1:numel(t)) = currentLevels;
% 
% 
% %===constant parameters===%
% %All of these can be found in Table 3
% gbar_K=36; gbar_Na=120; g_L=.3;
% E_K = -12; E_Na=115; E_L=10.6;
% C=1;
% 
% 
% %===set the initial states===%
% V=0; %Baseline voltage
% alpha_n = .01 * ( (10-V) / (exp((10-V)/10)-1) ); %Equation 12
% beta_n = .125*exp(-V/80); %Equation 13
% alpha_m = .1*( (25-V) / (exp((25-V)/10)-1) ); %Equation 20
% beta_m = 4*exp(-V/18); %Equation 21
% alpha_h = .07*exp(-V/20); %Equation 23
% beta_h = 1/(exp((30-V)/10)+1); %Equation 24
% 
% n(1) = alpha_n/(alpha_n+beta_n); %Equation 9
% m(1) = alpha_m/(alpha_m+beta_m); %Equation 18
% h(1) = alpha_h/(alpha_h+beta_h); %Equation 18
% 
% 
% for i=1:numel(t)-1 %Compute coefficients, currents, and derivates at each time step
%     
%     %---calculate the coefficients---%
%     %Equations here are same as above, just calculating at each time step
%     alpha_n(i) = .01 * ( (10-V(i)) / (exp((10-V(i))/10)-1) );
%     beta_n(i) = .125*exp(-V(i)/80);
%     alpha_m(i) = .1*( (25-V(i)) / (exp((25-V(i))/10)-1) );
%     beta_m(i) = 4*exp(-V(i)/18);
%     alpha_h(i) = .07*exp(-V(i)/20);
%     beta_h(i) = 1/(exp((30-V(i))/10)+1);
%     
%     
%     %---calculate the currents---%
%     I_Na = (m(i)^3) * gbar_Na * h(i) * (V(i)-E_Na); %Equations 3 and 14
%     I_K = (n(i)^4) * gbar_K * (V(i)-E_K); %Equations 4 and 6
%     I_L = g_L *(V(i)-E_L); %Equation 5
%     I_ion = I(i) - I_K - I_Na - I_L; 
%     
%     
%     %---calculate the derivatives using Euler first order approximation---%
%     V(i+1) = V(i) + deltaT*I_ion/C;
%     n(i+1) = n(i) + deltaT*(alpha_n(i) *(1-n(i)) - beta_n(i) * n(i)); %Equation 7
%     m(i+1) = m(i) + deltaT*(alpha_m(i) *(1-m(i)) - beta_m(i) * m(i)); %Equation 15
%     h(i+1) = h(i) + deltaT*(alpha_h(i) *(1-h(i)) - beta_h(i) * h(i)); %Equation 16
% 
% end
% 
% 
% V = V-70; %Set resting potential to -70mv
% 
% %===plot Voltage===%
% figure(6)
% for i = 1:4
%     subplot(2,2,i)
% plot(t,V,'LineWidth',3)
% hold on
% legend({'voltage'})
% ylabel('Voltage (mv)')
% xlabel('time (ms)')
% title('Voltage over Time in Simulated Neuron')
% end
% 
% 
% %===plot Conductance===%
% figure(7)
% p1 = plot(t,gbar_K*n.^4,'LineWidth',2);
% hold on
% p2 = plot(t,gbar_Na*(m.^3).*h,'r','LineWidth',2);
% legend([p1, p2], 'Conductance for Potassium', 'Conductance for Sodium')
% ylabel('Conductance')
% xlabel('time (ms)')
% title('Conductance for Potassium and Sodium Ions in Simulated Neuron')
% %Concentration of active ion channels multipled by a Ki that lies between a
% %range of 8-300 uM from paper of following html: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1171109/pdf/000131.pdf
% %on page 10/14. 
% % freq = zeros(1,7);
% 
figure(5)
for i = 1:4   
    subplot(2,2,i)
    %Iapp = (1 - outputHolder{i}(:,9))*(8*(1-outputHolder{i}(:,9)));
    Iapp = (outputHolder{i}(:,1) + outputHolder{i}(:,2))*80;
%     Iapp = 50;
    [t, yd] = actionPotential(timeCell{end}(end)*100, timeCell{end}, Iapp);
    freq(i) = freqCounter(yd, t, -50);
    plot(t, yd)
    xlabel('Time (ms)')
    ylabel('Membrane Voltage (mV)')
end

figure(6)
% Iapp = (outputHolder{5}(:,1) + outputHolder{5}(:,2))*80;
for i = 5:7   
    subplot(2,2,(i-4))
    Iapp = (outputHolder{i}(:,1) + outputHolder{i}(:,2))*(80);
    [t, yd] = actionPotential(timeCell{end}(end)*100, timeCell{end}, Iapp);
    freq(i) = freqCounter(yd, t, -50);
    plot(t, yd)
    xlabel('Time (ms)')
    ylabel('Membrane Voltage (mV)')
end

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