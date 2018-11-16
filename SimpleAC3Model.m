function SimpleAC3Model
cAMP(1) = .1*10^-6;
CNG(1) = .1*10^-6;
Ca(1) = .1*10^-6;
CaBP(1) = .1*10^-6;
CaCaM(1) = .1*10^-6;
dt = .01;

for t = 1:5000
    if 10 <= t && t <= 3000
        u = 200;
    else
        u = 0;
    end
    
    dcAMP = (2*.63*CNG(t)) - (2*.08*cAMP(t)^2*(.74-CNG(t))) - (3.16*cAMP(t)) - (47.02*cAMP(t)*CaCaM(t)) + u;
    dCNG = (.08*cAMP(t)^2*(.74-CNG(t))) - (.63*CNG(t)) - (163.17*CNG(t)*CaBP(t)^2);
    dCa = (47.29*CNG(t)) - (3.32*Ca(t)) - (.84*Ca(t)*(.74-CaBP(t))) + (.6*CaBP(t)) - (2*.01*Ca(t)^2*(1.3-CaCaM(t))) + (2*.1*CaCaM(t));
    dCaBP = (.84*Ca(t)*(.74-CaBP(t))) - (.6*CaBP(t));
    dCaCaM = (.01*Ca(t)^2*(1.3-CaCaM(t))) - (.1*CaCaM(t));
    
    cAMP(t+1) = cAMP(t) + dcAMP*dt;
    CNG(t+1) = CNG(t) + dCNG*dt;
    Ca(t+1) = Ca(t) + dCa*dt;
    CaBP(t+1) = CaBP(t) + dCaBP*dt;
    CaCaM(t+1) = CaCaM(t) + dCaCaM*dt;
    
    
end

subplot(5,1,1)
plot(0:.01:50,cAMP)
title('[cAMP] vs. t')

subplot(5,1,2)
plot(0:.01:50,CNG)
title('[CNG] vs. t')

subplot(5,1,3)
plot(0:.01:50,Ca)
title('[Ca] vs. t')

subplot(5,1,4)
plot(0:.01:50,CaBP)
title('[CaBP] vs. t')

subplot(5,1,5)
plot(0:.01:50,CaCaM)
title('[CaCaM] vs. t')
end