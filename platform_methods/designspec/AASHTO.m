function [Sa] = AASHTO(T,PGA,Ss,S1,siteclass)

switch siteclass
    case 'A'
        Fpga = 0.8;
        Fa   = 0.8;
        Fv   = 0.8;
    
    case 'B'
        Fpga = 1.0;
        Fa   = 1.0;
        Fv   = 1.0;
    
    case 'C'
        if     PGA<0.1,             Fpga = 1.2;
        elseif PGA>=0.1 && PGA<0.2, Fpga = 1.2+(1.2-1.2)/(0.2-0.1)*(PGA-0.1);
        elseif PGA>=0.2 && PGA<0.3, Fpga = 1.2+(1.1-1.2)/(0.3-0.2)*(PGA-0.2);
        elseif PGA>=0.3 && PGA<0.4, Fpga = 1.1+(1.0-1.1)/(0.4-0.3)*(PGA-0.3);
        elseif PGA>=0.4 && PGA<0.5, Fpga = 1.0+(1.0-1.0)/(0.5-0.4)*(PGA-0.4);
        elseif PGA>=0.5,            Fpga = 1.0;
        end
        if     Ss<0.25,             Fa = 1.2;
        elseif Ss>=0.25 && Ss<0.50, Fa = 1.2+(1.2-1.2)/(0.50-0.25)*(Ss-0.25);
        elseif Ss>=0.50 && Ss<0.75, Fa = 1.2+(1.1-1.2)/(0.75-0.50)*(Ss-0.50);
        elseif Ss>=0.75 && Ss<1.00, Fa = 1.1+(1.0-1.1)/(1.00-0.75)*(Ss-0.75);
        elseif Ss>=1.00 && Ss<1.25, Fa = 1.0+(1.0-1.0)/(1.25-1.00)*(Ss-1.00);
        elseif Ss>=1.25,            Fa = 1.0;
        end
        if     S1<0.1,              Fv = 1.7;
        elseif S1>=0.1 && S1<0.2,   Fv = 1.7+(1.6-1.7)/(0.2-0.1)*(S1-0.1);
        elseif S1>=0.2 && S1<0.3,   Fv = 1.6+(1.5-1.6)/(0.3-0.2)*(S1-0.2);
        elseif S1>=0.3 && S1<0.4,   Fv = 1.5+(1.4-1.5)/(0.4-0.3)*(S1-0.3);
        elseif S1>=0.4 && S1<0.5,   Fv = 1.4+(1.3-1.4)/(0.5-0.4)*(S1-0.4);
        elseif S1>=0.5,             Fv = 1.3;
        end
        
    case 'D'
        if     PGA<0.1,             Fpga = 1.6;
        elseif PGA>=0.1 && PGA<0.2, Fpga = 1.6+(1.4-1.6)/(0.2-0.1)*(PGA-0.1);
        elseif PGA>=0.2 && PGA<0.3, Fpga = 1.4+(1.2-1.4)/(0.3-0.2)*(PGA-0.2);
        elseif PGA>=0.3 && PGA<0.4, Fpga = 1.2+(1.1-1.2)/(0.4-0.3)*(PGA-0.3);
        elseif PGA>=0.4 && PGA<0.5, Fpga = 1.1+(1.0-1.1)/(0.5-0.4)*(PGA-0.4);
        elseif PGA>=0.5,            Fpga = 1.0;
        end
        if     Ss<0.25,             Fa = 1.6;
        elseif Ss>=0.25 && Ss<0.50, Fa = 1.6+(1.4-1.6)/(0.50-0.25)*(Ss-0.25);
        elseif Ss>=0.50 && Ss<0.75, Fa = 1.4+(1.2-1.4)/(0.75-0.50)*(Ss-0.50);
        elseif Ss>=0.75 && Ss<1.00, Fa = 1.2+(1.1-1.2)/(1.00-0.75)*(Ss-0.75);
        elseif Ss>=1.00 && Ss<1.25, Fa = 1.1+(1.0-1.1)/(1.25-1.00)*(Ss-1.00);
        elseif Ss>=1.25,            Fa = 1.0;
        end
        if     S1<0.1,              Fv = 2.4;
        elseif S1>=0.1 && S1<0.2,   Fv = 2.4+(2.0-2.4)/(0.2-0.1)*(S1-0.1);
        elseif S1>=0.2 && S1<0.3,   Fv = 2.0+(1.8-2.0)/(0.3-0.2)*(S1-0.2);
        elseif S1>=0.3 && S1<0.4,   Fv = 1.8+(1.6-1.8)/(0.4-0.3)*(S1-0.3);
        elseif S1>=0.4 && S1<0.5,   Fv = 1.6+(1.5-1.6)/(0.5-0.4)*(S1-0.4);
        elseif S1>=0.5,             Fv = 1.5;
        end
        
    case 'E'
        if     PGA<0.1,             Fpga = 2.5;
        elseif PGA>=0.1 && PGA<0.2, Fpga = 2.5+(1.7-2.5)/(0.2-0.1)*(PGA-0.1);
        elseif PGA>=0.2 && PGA<0.3, Fpga = 1.7+(1.2-1.7)/(0.3-0.2)*(PGA-0.2);
        elseif PGA>=0.3 && PGA<0.4, Fpga = 1.2+(0.9-1.2)/(0.4-0.3)*(PGA-0.3);
        elseif PGA>=0.4 && PGA<0.5, Fpga = 0.9+(0.9-0.9)/(0.5-0.4)*(PGA-0.4);
        elseif PGA>=0.5,            Fpga = 0.9;
        end
        if     Ss<0.25,             Fa = 2.5;
        elseif Ss>=0.25 && Ss<0.50, Fa = 2.5+(1.7-2.5)/(0.50-0.25)*(Ss-0.25);
        elseif Ss>=0.50 && Ss<0.75, Fa = 1.7+(1.2-1.7)/(0.75-0.50)*(Ss-0.50);
        elseif Ss>=0.75 && Ss<1.00, Fa = 1.2+(0.9-1.2)/(1.00-0.75)*(Ss-0.75);
        elseif Ss>=1.00 && Ss<1.25, Fa = 0.9+(0.9-0.9)/(1.25-1.00)*(Ss-1.00);
        elseif Ss>=1.25,            Fa = 0.9;
        end  
        if     S1<0.1,              Fv = 3.5;
        elseif S1>=0.1 && S1<0.2,   Fv = 3.5+(3.2-3.5)/(0.2-0.1)*(S1-0.1);
        elseif S1>=0.2 && S1<0.3,   Fv = 3.2+(2.8-3.2)/(0.3-0.2)*(S1-0.2);
        elseif S1>=0.3 && S1<0.4,   Fv = 2.8+(2.4-2.8)/(0.4-0.3)*(S1-0.3);
        elseif S1>=0.4 && S1<0.5,   Fv = 2.4+(2.4-2.4)/(0.5-0.4)*(S1-0.4);
        elseif S1>=0.5,             Fv = 2.4;
        end    
end

As  = Fpga*PGA;
Sds = Fa*Ss;
Sd1 = Fv*S1;
Ts  = Sd1/Sds;
To  = 0.2*Ts;

np = length(T);
Sa = ones(1,np);

for i=1:np
    if     T(i)<=To,             Sa(i)=As+(Sds-As)*(T(i)/To);
    elseif To<=T(i) && T(i)<=Ts, Sa(i)=Sds;
    elseif T(i)>Ts,              Sa(i)=Sd1/T(i);
    end
end

end

