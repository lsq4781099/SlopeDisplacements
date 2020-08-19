function [Sa] = NTCD2004(T,zone)

switch zone
    case 'I',    c = 0.16; ao = 0.04; Ta = 0.20; Tb = 1.35; r = 1.00;
    case 'II',   c = 0.32; ao = 0.08; Ta = 0.20; Tb = 1.35; r = 1.33;
    case 'IIIa', c = 0.40; ao = 0.10; Ta = 0.53; Tb = 1.80; r = 2.00;
    case 'IIIb', c = 0.45; ao = 0.11; Ta = 0.85; Tb = 3.00; r = 2.00;
    case 'IIIc', c = 0.40; ao = 0.10; Ta = 1.25; Tb = 4.20; r = 2.00;
    case 'IIId', c = 0.30; ao = 0.10; Ta = 0.85; Tb = 4.20; r = 2.00;
end

np = length(T);
Sa = ones(1,np)*ao;

for i=1:np
    if     Ta>=T(i)           , Sa(i)=ao+(c-ao)*T(i)/Ta;
    elseif Ta<=T(i) && T(i)<Tb, Sa(i)=c;
    elseif Tb<=T(i)           , Sa(i)=c*(Tb/T(i))^r;
    end
end

end

