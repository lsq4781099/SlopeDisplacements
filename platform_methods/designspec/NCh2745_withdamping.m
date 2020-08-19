function [Sa] = NCh2745(T,zone,soiltype,beta)

switch zone
    case 1, Z = 3/4;
    case 2, Z = 4/4;
    case 3, Z = 5/4;
end

switch soiltype
    case 'I',   A = 0.40; Ta = 0.03; Tb = 0.110; Tc = 0.29; Td = 2.51; aA = 1085; aV = 50;  aD = 20;
    case 'II',  A = 0.41; Ta = 0.03; Tb = 0.200; Tc = 0.54; Td = 2.00; aA = 1100; aV = 94;  aD = 30;
    case 'III', A = 0.45; Ta = 0.03; Tb = 0.375; Tc = 0.68; Td = 1.58; aA = 1212; aV = 131; aD = 33;
end

if     beta<=0.02,                B = 0.65;
elseif beta> 0.02  && beta<=0.05, B = 0.65+(1.00-0.65)/(0.05-0.02)*(beta-0.02);
elseif beta> 0.05  && beta<=0.10, B = 1.00+(1.37-1.00)/(0.10-0.05)*(beta-0.05);
elseif beta> 0.10  && beta<=0.15, B = 1.37+(1.67-1.37)/(0.15-0.10)*(beta-0.10);
elseif beta> 0.15  && beta<=0.20, B = 1.67+(1.94-1.67)/(0.20-0.15)*(beta-0.15);
elseif beta> 0.20  && beta<=0.25, B = 1.94+(2.17-1.94)/(0.25-0.20)*(beta-0.20);
elseif beta> 0.25  && beta<=0.30, B = 2.17+(2.38-2.17)/(0.30-0.25)*(beta-0.25);
elseif beta> 0.30  && beta<=0.50, B = 2.38+(3.02-2.38)/(0.50-0.30)*(beta-0.30);
elseif beta> 0.50,                B = 3.02;
end

A  = A*981;
np = length(T);
Sa = ones(1,np)*A;

for i=1:np
    if     Ta<=T(i) && T(i)<Tb, Sa(i)=A+(aA-A)/(Tb-Ta)*(T(i)-Ta);
    elseif Tb<=T(i) && T(i)<Tc, Sa(i)=aA;
    elseif Tc<=T(i) && T(i)<Td, Sa(i)=2*pi/T(i)*aV;
    elseif Td<=T(i)           , Sa(i)=4*pi^2/T(i)^2*aD;
    end
end

Sa = Sa/981*Z/B;

end

