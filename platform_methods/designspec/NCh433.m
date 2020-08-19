function [Sa] = NCh433(T,category,zone,soiltype)

switch category
    case 'A', I = 1.2;
    case 'B', I = 1.2;
    case 'C', I = 1.0;
    case 'D', I = 0.6;
end

switch zone
    case 1, Ao = 0.2;
    case 2, Ao = 0.3;
    case 3, Ao = 0.4;
end

switch soiltype
    case 'I'  , To = 0.15; p = 2.0;
    case 'II' , To = 0.30; p = 1.5;
    case 'III', To = 0.75; p = 1.0;
    case 'IV' , To = 1.20; p = 1.0;
end

a  = (1+4.5*(T./To).^p)./(1+(T./To).^3);
Sa = I*Ao*a;

end

