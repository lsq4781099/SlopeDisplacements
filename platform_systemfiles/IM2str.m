function[str]=IM2str(T)

str      = cell(size(T));
SPECTRAL = imag(T)~=0;

for i=1:numel(T)
    val = T(i);
    
    if SPECTRAL(i)==0
        if val==0,      str{i} = 'PGA';
        elseif val==-1, str{i} = 'PGV';
        elseif val==-2, str{i} = 'PGD';
        elseif val==-3, str{i} = 'Duration';
        elseif val==-4, str{i} = 'CAV';
        elseif val==-5, str{i} = 'AI';
        elseif val==-6, str{i} = 'VGI';
        elseif val>0  , str{i} = sprintf('SA(T=%g)',val);
        end
    end
    
    if SPECTRAL(i)==1
        RPART = real(val);
        switch imag(val)
            case 1, str{i} = sprintf('SV(T=%g)',RPART);
            case 2, str{i} = sprintf('SD(T=%g)',RPART);
            case 3, str{i} = sprintf('H/V(T=%g)',RPART);
        end
    end
    
end
