function [Sa] = SENCICO(T,category,zone,soiltype)

switch category
    case 'A1', U = 1.5; 
        if zone==3 || zone==4, U = 1.0; end
    case 'A2', U = 1.5;
    case 'B',  U = 1.3;
    case 'C',  U = 1.0;
    case 'D',  U = 1.0;
end

switch zone
    case 1, Z = 0.10;
    case 2, Z = 0.25;
    case 3, Z = 0.35;
    case 4, Z = 0.45;
end

switch soiltype
    case 'S0', Tp = 0.30; Tl = 3.00;
    case 'S1', Tp = 0.40; Tl = 2.50;
    case 'S2', Tp = 0.60; Tl = 2.00;
    case 'S3', Tp = 1.00; Tl = 1.60;
end

if     Z==0.10 && Tp==0.3, S = 0.80;
elseif Z==0.10 && Tp==0.4, S = 1.00;
elseif Z==0.10 && Tp==0.6, S = 1.60;
elseif Z==0.10 && Tp==1.0, S = 2.00;
elseif Z==0.25 && Tp==0.3, S = 0.80;
elseif Z==0.25 && Tp==0.4, S = 1.00;
elseif Z==0.25 && Tp==0.6, S = 1.20;
elseif Z==0.25 && Tp==1.0, S = 1.40;
elseif Z==0.35 && Tp==0.3, S = 0.80;
elseif Z==0.35 && Tp==0.4, S = 1.00;
elseif Z==0.35 && Tp==0.6, S = 1.15;
elseif Z==0.35 && Tp==1.0, S = 1.20;
elseif Z==0.45 && Tp==0.3, S = 0.80;
elseif Z==0.45 && Tp==0.4, S = 1.00;
elseif Z==0.45 && Tp==0.6, S = 1.05;
elseif Z==0.45 && Tp==1.0, S = 1.10;
end

np = length(T);
C  = ones(1,np);

for i=1:np
    if     T(i)<=Tp,           C(i) = 2.5;
    elseif Tp<T(i) && T(i)<Tl, C(i) = 2.5*(Tp/T(i));
    elseif T(i)>=Tl,           C(i) = 2.5*(Tp*Tl/T(i)^2);
    end
end

g  = 1;
Sa = (Z*U*C*S)*g;

end

