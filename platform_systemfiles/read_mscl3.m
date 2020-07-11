function[mscl]=read_mscl3(t)

mscl=createObj('mscl');
if isempty(t)
    return
end

nt=length(t);
source = cell(nt,1);
adp    = cell(nt,1);
num    = zeros(nt,4);
strpat ='%s %s %s %f %f';

for i=1:nt
    D=textscan(t{i},strpat);
    source{i}=D{2}{1};
    adp{i} = D{3}{1};
    num(i,:)=[1,nan,horzcat(D{4:5})]; % the "1" means NMmin
end

mscl.source=source;
mscl.adp=adp;
mscl.num=num;