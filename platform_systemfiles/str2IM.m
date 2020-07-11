function[IM]=str2IM(str)

str = upper(str);
if ~iscell(str)
    str = {str};
end

NIM = length(str);
IM  = nan(size(str));
SPECTRAL = contains(str,{'SA','SV','SD','H/V'});

for i=1:NIM
    if SPECTRAL(i)==0
        switch str{i}
            case 'PGA',  IM(i)=0;
            case 'PGV',  IM(i)=-1;
            case 'PGD',  IM(i)=-2;
            case 'DURATION',IM(i)=-3;
            case 'CAV',  IM(i)=-4;
            case 'AI' ,  IM(i)=-5;
            case 'VGI',  IM(i)=-6;
        end
    end
    
    if SPECTRAL(i)==1
        st = str{i};
        st = strrep(st,'TS','');
        st = strrep(st,'(T','');
        st = strrep(st,')','');
        st = regexp(st,'\=','split');
        switch st{1}
            case 'SA'  , IM(i)=str2double(st{2});
            case 'SV'  , IM(i)=str2double(st{2})+1i;
            case 'SD'  , IM(i)=str2double(st{2})+2i;
            case 'H/V' , IM(i)=str2double(st{2})+3i;
        end
    end
end