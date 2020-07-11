function [NMmin,param,val] = mMRgetparam(handles) %#ok<*INUSD,*DEFNU>

ch=get(handles.panel2,'children');
ch2=ch(handles.edit);
val=vertcat(ch2.Value);
str=cell(length(ch2),1);
for i=1:length(ch2)
    str{i}=ch2(i).String;
    if val(i)~=0 && size(str{i},1)>1
        str{i}=str{i}{val(i)};
    end
end
val = [val;handles.MRselect.Value];

methods = pshatoolbox_methods(2);

switch methods(handles.MRselect.Value).str
    case 'delta'
        NMmin    = str2double(str{1});
        M        = str2double(str{2});
        param    = M;
        
    case 'truncexp'
        NMmin    = str2double(str{1});
        bvalue   = str2double(str{2});
        Mmin     = str2double(str{3});
        Mmax     = str2double(str{4});
        param    = [bvalue,Mmin,Mmax];
        
    case 'truncnorm'
        NMmin    = str2double(str{1});
        Mmin     = str2double(str{2});
        Mmax     = str2double(str{3});
        Mchar    = str2double(str{4});
        sigmaM   = str2double(str{5});
        param    = [Mmin,Mmax,Mchar,sigmaM];
        
    case 'yc1985'
        NMmin    = str2double(str{1});
        bvalue   = str2double(str{2});
        Mmin     = str2double(str{3});
        Mchar    = str2double(str{4});
        param    = [bvalue,Mmin,Mchar];
        
    case 'magtable'
        Mmin     = str2double(str{1});
        binwidth = str2double(str{2});
        occurrance = ch2(14).String;
        occurrance = str2double(occurrance);
        NMmin    = occurrance(1);
        param    = [Mmin,binwidth,occurrance(:)'];
end
