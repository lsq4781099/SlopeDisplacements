function [Sa] = MdC(T,zone,soiltype,CI)

switch zone
    case 1, Ao = 0.2;
    case 2, Ao = 0.3;
    case 3, Ao = 0.4;
end

switch soiltype
    case 'I'  , S = 0.9; T1 = 0.20; K2 = 0.513;
    case 'II' , S = 1.0; T1 = 0.30; K2 = 0.672;
    case 'III', S = 1.2; T1 = 0.70; K2 = 1.182;
    case 'IV' , S = 1.3; T1 = 1.10; K2 = 1.598;
end

switch CI
    case 'I' , K1 = 1.0;
    case 'II', K1 = 0.8;
end

np = length(T);
Sa = ones(1,np);

for i=1:np
    if      T(i)<= T1,   Sa(i) = 1.5*K1*S*Ao;
    elseif  T1  <  T(i), Sa(i) = K1*K2*S*Ao/T(i)^(2/3);
    end
end

end

