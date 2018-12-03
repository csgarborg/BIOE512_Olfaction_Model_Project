function [timeOutput, APOutput] = actionPotential(simulationTime, orgTime, current)
 
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
    
    V=0; %Baseline voltage
    alpha_n = .01 * ( (10-V) / (exp((10-V)/10)-1) ); %Potassium Channel gate
    beta_n = .125*exp(-V/80); %Potassium Channel gate
    alpha_m = .1*( (25-V) / (exp((25-V)/10)-1) ); %Sodium channel gate
    beta_m = 4*exp(-V/18); %Sodium channel gate
    alpha_h = .07*exp(-V/20); %Sodium channel gate (Ball and socket)
    beta_h = 1/(exp((30-V)/10)+1); %Sodium channel gate (Ball and socket)

    %Initializing the probability values of each ion gate being opened at start     
    n(1) = alpha_n/(alpha_n+beta_n);
    m(1) = alpha_m/(alpha_m+beta_m);
    h(1) = alpha_h/(alpha_h+beta_h); 
    
    for i=1:numel(t)-1 %Compute coefficients, currents, and derivates at each time step
        alpha_n(i) = .01 * ( (10-V(i)) / (exp((10-V(i))/10)-1) );
        beta_n(i) = .125*exp(-V(i)/80);
        alpha_m(i) = .1*( (25-V(i)) / (exp((25-V(i))/10)-1) );
        beta_m(i) = 4*exp(-V(i)/18);
        alpha_h(i) = .07*exp(-V(i)/20);
        beta_h(i) = 1/(exp((30-V(i))/10)+1);
        
        %calculate the currents
        I(i) = interp1(orgTime, current, t(i)/100);
        I_Na = (m(i)^3) * gbar_Na * h(i) * (V(i)-E_Na); 
        I_K = (n(i)^4) * gbar_K * (V(i)-E_K); 
        I_L = g_L *(V(i)-E_L); 
        I_ion = I(i) - I_K - I_Na - I_L;
        
        V(i+1) = V(i) + deltaT*I_ion/C;
        n(i+1) = n(i) + deltaT*(alpha_n(i) *(1-n(i)) - beta_n(i) * n(i));
        m(i+1) = m(i) + deltaT*(alpha_m(i) *(1-m(i)) - beta_m(i) * m(i));
        h(i+1) = h(i) + deltaT*(alpha_h(i) *(1-h(i)) - beta_h(i) * h(i));     
    end
    
    APOutput = V-70; %Setting all values with respect to resting voltage
    timeOutput = t;
end
