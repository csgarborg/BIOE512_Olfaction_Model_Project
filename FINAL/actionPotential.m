function [timeOutput, APOutput, gbar_K, gbar_Na, m, h, n] = actionPotential(simulationTime, current)
%action potential model based off the Hodgkin and Huxley model for action
%potential stimulation in Squid. This code is also based off the code from
%Andrew Jahn. The full reference is listed below:
%Jahn, Andrew. “Introduction to Computational Modeling: Hodgkin-Huxley Model.
%” Andy's Brain Blog, 15 Oct. 2013, andysbrainblog.blogspot.com/2013/10/...
% ... introduction-to-computational-modeling.html.

    %Initializing variables and resting voltages for ions
    gbar_K=36; 
    gbar_Na=120; 
    g_L=.3;
    E_K = -12;
    E_Na=115;
    E_L=10.6;
    C=1;
    
    deltaT = 0.01;
    t=0:deltaT:simulationTime;
    V = zeros(1,length(t));
    n = zeros(1,length(t));
    m = zeros(1,length(t));
    h = zeros(1,length(t));
    alpha_n = zeros(1,length(t));
    beta_n = zeros(1,length(t));
    alpha_m = zeros(1,length(t));
    beta_m = zeros(1,length(t));
    alpha_h = zeros(1,length(t));
    beta_h = zeros(1,length(t));
    I = zeros(1, length(t));
    I_ion = I;
    I_Na = I;
    I_K = I;
    I_L = I;
    
    V(1)=0; %Baseline voltage
    alpha_n(1) = .01 * ( (10-V(1)) / (exp((10-V(1))/10)-1) ); %Potassium Channel gate
    beta_n(1) = .125*exp(-V(1)/80); %Potassium Channel gate
    alpha_m(1) = .1*( (25-V(1)) / (exp((25-V(1))/10)-1) ); %Sodium channel gate
    beta_m(1) = 4*exp(-V(1)/18); %Sodium channel gate
    alpha_h(1) = .07*exp(-V(1)/20); %Sodium channel gate (Ball and socket)
    beta_h(1) = 1/(exp((30-V(1))/10)+1); %Sodium channel gate (Ball and socket)

    %Initializing the probability values of each ion gate being opened at start     
    n(1) = alpha_n(1)/(alpha_n(1)+beta_n(1));
    m(1) = alpha_m(1)/(alpha_m(1)+beta_m(1));
    h(1) = alpha_h(1)/(alpha_h(1)+beta_h(1)); 
    
    %Compute coefficients, currents, and derivates at each time step
    for i=1:numel(t)-1 
        alpha_n(i) = .01 * ( (10-V(i)) / (exp((10-V(i))/10)-1) );
        beta_n(i) = .125*exp(-V(i)/80);
        alpha_m(i) = .1*( (25-V(i)) / (exp((25-V(i))/10)-1) );
        beta_m(i) = 4*exp(-V(i)/18);
        alpha_h(i) = .07*exp(-V(i)/20);
        beta_h(i) = 1/(exp((30-V(i))/10)+1);
        
        %calculate the currents
        I(i) = current(i);
        I_Na(i) = (m(i)^3) * gbar_Na * h(i) * (V(i)-E_Na); 
        I_K(i) = (n(i)^4) * gbar_K * (V(i)-E_K); 
        I_L(i) = g_L *(V(i)-E_L); 
        I_ion(i) = I(i) - I_K(i) - I_Na(i) - I_L(i);
        
        V(i+1) = V(i) + deltaT*I_ion(i)/C;
        n(i+1) = n(i) + deltaT*(alpha_n(i) *(1-n(i)) - beta_n(i) * n(i));
        m(i+1) = m(i) + deltaT*(alpha_m(i) *(1-m(i)) - beta_m(i) * m(i));
        h(i+1) = h(i) + deltaT*(alpha_h(i) *(1-h(i)) - beta_h(i) * h(i));     
    end
    
    APOutput = V-70; %Setting all values with respect to resting voltage
    timeOutput = t;
    fprintf('done \n')
end
