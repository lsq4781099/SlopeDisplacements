function[str]=field2str(str)

% fixes strings of IMs
N = length(str);
str = strrep(str,' ','');
NOTNAN=~isnan(str2double(str));
for i=1:N
    if NOTNAN(i)
        switch str{i}
            case '0'   , str{i}='PGA';
            case '-1'  , str{i}='PGV';
            case '-2'  , str{i}='PGD';
            case '-3'  , str{i}='Duration';
            case '-4'  , str{i}='CAV';
            case '-5'  , str{i}='AI';
            case '-6'  , str{i}='VGI';
            otherwise
                if     contains(str{i},'+1i'),  num = real(str2double(str{i})); str{i}=sprintf('SV(T=%g)' ,num);
                elseif contains(str{i},'+2i'),  num = real(str2double(str{i})); str{i}=sprintf('SD(T=%g)' ,num);
                elseif contains(str{i},'+3i'),  num = real(str2double(str{i})); str{i}=sprintf('H/V(T=%g)',num);
                else
                    num = real(str2double(str{i}));
                    if num>0
                        str{i}=sprintf('SA(T=%g)' ,num);
                    else
                        str{i}='Error';
                    end
                end
        end
    end
end
