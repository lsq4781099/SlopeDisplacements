function [Sa] = NCh2369(T,category,zone,soiltype,xi,R)

switch category
    case 'C1', I = 1.2;
    case 'C2', I = 1.0;
    case 'C3', I = 0.8;
end

switch zone
    case 1, Ao = 0.2; factor = 0.50;
    case 2, Ao = 0.3; factor = 0.75;
    case 3, Ao = 0.4; factor = 1.00;
end

switch soiltype
    case 'I'  , Tp = 0.20; n = 1.00;
    case 'II' , Tp = 0.35; n = 1.33;
    case 'III', Tp = 0.62; n = 1.80;
    case 'IV' , Tp = 1.35; n = 1.80;
end

if     R==1 && xi==0.02, Cmax = 0.79;
elseif R==1 && xi==0.03, Cmax = 0.68;
elseif R==1 && xi==0.05, Cmax = 0.55;
elseif R==2 && xi==0.02, Cmax = 0.60;
elseif R==2 && xi==0.03, Cmax = 0.49;
elseif R==2 && xi==0.05, Cmax = 0.42;
elseif R==3 && xi==0.02, Cmax = 0.40;
elseif R==3 && xi==0.03, Cmax = 0.34;
elseif R==3 && xi==0.05, Cmax = 0.28;
elseif R==4 && xi==0.02, Cmax = 0.32;
elseif R==4 && xi==0.03, Cmax = 0.27;
elseif R==4 && xi==0.05, Cmax = 0.22;
elseif R==5 && xi==0.02, Cmax = 0.26;
elseif R==5 && xi==0.03, Cmax = 0.23;
elseif R==5 && xi==0.05, Cmax = 0.18;
end

Cmax = Cmax*factor;
g    = 1;
cap  = I*Cmax*g;
Sa   = min(2.75*Ao*I*(Tp./T).^n*(0.05/xi)^0.4,cap);

end

