function UN=addIMunits(str)

if isnumeric(str)
    str = IM2str(str);
    str = str{1};
end

UN = str; % this takes care of dimensionless IMs
switch str
    case 'PGA'         ,UN='PGA (g)'; return
    case 'PGV'         ,UN='PGV (cm/s)';return
    case 'PGD'         ,UN='PGD (cm)';return
    case 'Duration'    ,UN='Duration (s)';return
    case 'CAV'         ,UN='CAV (cm*s)';return
    case 'AI'          ,UN='AI (cm/s)';return
    case 'VGI'         ,UN='VGI (cm/s)';return
end

if contains(str,'SA('), UN = [str, ' (g)']; return; end
if contains(str,'SV('), UN = [str, ' (cm/s)']; return; end
if contains(str,'SD('), UN = [str, ' (cm)']; return; end
if contains(str,'SD('), UN = [str, ' (cm)']; return; end

