function[param]=loadsiteLIBS(filename)

fid  = fopen(filename);
if fid==-1
    return
end
data = textscan(fid,'%s','delimiter','\n');
data = data{1};
fclose(fid);

ind=strfind(data,'#');
for i=1:size(data,1)
    if ~isempty(ind{i})
        II = ind{i}(1);
        data{i}(II:end)=[];
    end
    data{i}=deblank(data{i});
end

% removes empty lines
emptylist= [];
for i=1:size(data,1)
    if isempty(data{i,1})
        emptylist=[emptylist,i]; %#ok<*AGROW>
    end
end
data(emptylist,:)=[];

% removes multiple spaces
data=regexprep(data,' +',' ');

B     = regexp(data{1} ,'\:','split'); param.B     = str2double(B{2});
L     = regexp(data{2} ,'\:','split'); param.L     = str2double(L{2});
Df    = regexp(data{3} ,'\:','split'); param.Df    = str2double(Df{2});
Q     = regexp(data{4} ,'\:','split'); param.Q     = eval(['[',Q{2},']']);
type  = regexp(data{5} ,'\:','split'); param.type  = strtrim(type{2});
wt    = regexp(data{6} ,'\:','split'); param.wt    = str2double(wt{2});
LPC   = regexp(data{7} ,'\:','split'); param.LPC   = strtrim(LPC{2});
th1   = regexp(data{8} ,'\:','split'); param.th1   = str2double(th1{2});
th2   = regexp(data{9} ,'\:','split'); param.th2   = str2double(th2{2});
N1    = regexp(data{10},'\:','split'); param.N1    = str2double(N1{2});
N2    = regexp(data{11},'\:','split'); param.N2    = str2double(N2{2});
meth  = regexp(data{12},'\:','split'); param.meth  = strtrim(meth{2});
N160  = regexp(data{13},'\:','split'); param.N160  = str2double(regexp(strtrim( N160{2}),'\ ','split'));
qc1N  = regexp(data{14},'\:','split'); param.qc1N  = str2double(regexp(strtrim( qc1N{2}),'\ ','split'));
thick = regexp(data{15},'\:','split'); param.thick = str2double(regexp(strtrim(thick{2}),'\ ','split'));
d2mat = regexp(data{16},'\:','split'); param.d2mat = str2double(regexp(strtrim(d2mat{2}),'\ ','split'));

param.CPT = regexp(data(18:end,:),'\ ','split');
param.CPT = str2double(vertcat(param.CPT{:}));
