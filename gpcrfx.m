function yd = gpcrfx(~,y)

    global a2 a3 a4 a5 a6 a7 a8;

    %%Initial Conditions 
    Tri_GP = y(1);
    B_gamma = y(2);
    Ga_GTP = y(3);
    Ga_GDP = y(4);

    %%Constants (global variables)
    GAPs = a2;
    Kass = a3;
    Kdiss = a4;
    k2 = a5;
    RC = a6;
    Khydro = a7;
    k3 = a8;
    
    %ODE system for G protein system
    v1 = Kass*Ga_GDP*B_gamma;
    v2 = (Kdiss*Tri_GP*RC)/(k2+Tri_GP);
    v3 = (Khydro*Ga_GTP*GAPs)/(k3+Ga_GTP);

    yd(1) = v1 - v2;
    yd(2) = v2 - v1;
    yd(3) = v2 - v3;
    yd(4) = v3 - v1;
    yd = yd';
    
end

