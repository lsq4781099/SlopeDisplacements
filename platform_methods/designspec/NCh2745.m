function [Sa] = NCh2745(T,zone,soiltype)

switch zone
    case 1, Z=3/4;
    case 2, Z=4/4;
    case 3, Z=5/4;
end

switch soiltype
    case 'I'  , A = 0.40; Ta = 0.03; Tb = 0.110; Tc = 0.29; Td = 2.51; aA = 1085; aV = 50;  aD = 20;
    case 'II' , A = 0.41; Ta = 0.03; Tb = 0.200; Tc = 0.54; Td = 2.00; aA = 1100; aV = 94;  aD = 30;
    case 'III', A = 0.45; Ta = 0.03; Tb = 0.375; Tc = 0.68; Td = 1.58; aA = 1212; aV = 131; aD = 33;
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

Sa = Z*Sa/981;

end

