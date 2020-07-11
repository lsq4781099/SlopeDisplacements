function[mscl]=read_mscl2(t)

mscl=createObj('mscl');
if isempty(t)
    return
end

nt = size(t,1);
source = cell(nt,1);
num    = cell(nt,1);

for i=1:nt
    ti        = regexp(t{i},'\ ','split');
    source{i} = ti{2};
    num{i}    = [1,str2double(ti{5}),str2double(ti(3:end))];
end

mscl.source=source;
mscl.num=num;
