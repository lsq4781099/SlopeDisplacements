function mGMPEcheck_param(handles)

val   = handles.GMPEselect.Value;
isSUB = any(handles.methods(val).mech(1:2));

ch  = get(handles.panel2,'children');
txt = ch(handles.text);
edi = ch(handles.edit);

%% makes sure numbers are positive
for i=1:length(txt)
    switch txt(i).String
        case 'M'   , edi(i).String=strrep(edi(i).String,'-','');
        case 'Rrup', edi(i).String=strrep(edi(i).String,'-','');
        case 'Rhyp', edi(i).String=strrep(edi(i).String,'-','');
        case 'Rjb' , edi(i).String=strrep(edi(i).String,'-','');
        case 'Rx'  
        case 'Ry0' , edi(i).String=strrep(edi(i).String,'-','');
        case 'Zhyp', edi(i).String=strrep(edi(i).String,'-','');
        case 'Ztor', edi(i).String=strrep(edi(i).String,'-','');
        case 'Zbot', edi(i).String=strrep(edi(i).String,'-','');
        case 'Vs30', edi(i).String=strrep(edi(i).String,'-','');
        case 'Z10' , edi(i).String=strrep(edi(i).String,'-','');
        case 'Z25' , edi(i).String=strrep(edi(i).String,'-','');
    end
end

%% Check that rrup>ztor or zhyp for subduction GMMs

if isSUB
    str = {txt.String};
    [~,ind1]=intersect(str,'Rrup');
    [~,ind2]=intersect(str,'Zhyp');
    [~,ind3]=intersect(str,'Ztor');
    
    if ~isempty(ind2)
        Rrup = str2double(edi(ind1).String);
        Zhyp = str2double(edi(ind2).String);
        if Rrup<Zhyp
            edi(ind1).String=edi(ind2).String;
        end
        return
    end
    
    
    if ~isempty(ind3)
        Rrup = str2double(edi(ind1).String);
        Ztor = str2double(edi(ind3).String);
        if Rrup<Ztor
            edi(ind1).String=edi(ind3).String;
        end
        return
    end
end