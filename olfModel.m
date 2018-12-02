function output = olfModel(time,input)

    global od rc os os2 pka AKmax ak2 PKmax pk2 go rg rg2 ac a ckk A3max A2max...
    ca apd PDmax cng CNmax cn2 cam cc cc2 CMmax ifa ef cac CLmax...
    cl2 cpd CPmax cp2 CDmax CKmax ck2 K1 K2 K3 N1 N2 N3 tempHold;
    
    %Setting concentrations
    x = input(1); %concentration of active CNG channels
    y = input(2); %concentration of active Ca-Cl channels
    z = input(3); %concentration of odorant-receptor complex
    w = input(4); %concentration of active Golf
    v = input(5); %concentration of active AC
    u = input(6); %concentration of cAMO 
    s = input(7); %concentration of active PKA
    r = input(8); %concentration of inactive odorant-receptor complex
    q = input(9); %concentration of intracellular Ca
    p = input(10); %concentration of Ca4CAM
    o = input(11); %concentration of CAM-PDE
    n = input(12); %concentration of active CAM-kinase II
    
    %Setting duration of odorant exposure (starts at time = 0), controls
    %pulse width of odorant exposure and number of exposures 
    duration = 10;%seconds
    numberOfRepeats = 1;
    delay = 0; %seconds
    for i = 1:numberOfRepeats
       if time <= duration
           od = tempHold; %constant from global variable
           break;
       elseif (duration*i + delay) <= time &&  time <= (duration*i+2*delay) 
           od = tempHold;
           break;
       else
           od = 0;
       end
    end
    
    %Initializing Variables
    io = cpd-o; %concentration of inactive CAM-PDE
    in = ckk-n; %concentration of inactive CAM-KINASE II
    iv = ac-v; %concentration of inactive AC
    ix = cng-x; %concentration of inacive CNG channels
    iy = cac-y; %concentration of inactive Ca-Cl channels
    iw = go-w; % concentration of inactive Golf
    urc = rc-z-r; %concentration of unbound receptors 
    uod = abs(od-z-r); %concentration of unbound oderants 
    is = pka-s; %concentration of inactive PKA
    cd = CDmax*o; %rate of cAMP hydrolysis by CAM-PDE
    cp = CPmax*p; %rate of CAM-PDE activation by Ca4CAM
    ck = CKmax*p; %rate of CAM-kinase II activation by Ca4CAM
    cm = CMmax*p; %rate of CNG channel inactivation by Ca4CAM
    a2 = A2max*n; % rate of AC inactivation by CAM-KINASE II
    pk = PKmax*s; % rate of receptor inactivation by active PKA
    ak = AKmax*u; %rate of PKA activation by cAMP
    pd = PDmax*apd; % rate of cAMP hydrolysis by cAMP-PDE
    cn = CNmax*u^N1/(u^N1+K1^N1); %rate of CNG channel activation by cAMP
    cl = CLmax*q^N2/(q^N2+K2^N2); %rate of Ca-Cl channel activation by Ca ions
    a3 = A3max*q^N3/(q^N3+K3^N3); %rate of AC inactivation by Ca ions
    
    %ODE Systems
    output(1) = cn*ix-cn2*x-cm*x; %dx/dt
    output(2) = cl*iy-cl2*y; %dy/dt
    output(3) = urc*os*uod-os2*z-pk*z+pk2*r; %dz/dt
    output(4) = rg*iw*z-rg2*w; %dw/dt
    output(5) = iv*a*w-a3*v-a2*v; %dv/dt
    output(6) = ca*v-pd*u-cd*u; %du/dt
    output(7) = is*ak-ak2*s; %ds/dt
    output(8) = pk*z-pk2*r; %dr/dt
    output(9) = ifa*x-ef*q; %dq/dt
    output(10) = cc*cam*(q/4)-cc2*p; %dp/dt
    output(11) = cp*io-cp2*o; %do/dt
    output(12) = ck*in-ck2*n; %dn/dt
    
    output = output';
end
