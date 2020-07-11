function[mscl]=read_mscl1(t,mtype)

mscl=createObj('mscl');
if isempty(t)
    return
end

st = sprintf('%s ',mtype);
t  = strrep(t,st,'');
t  = strrep(t,' NM ',' 1 ');
t  = strrep(t,' SR ',' 2 ');

nt=length(t);
source = cell(nt,1);

switch mtype
    case 'delta',     num = zeros(nt,3); strpat='%s %f %f %f';
    case 'truncexp',  num = zeros(nt,5); strpat='%s %f %f %f %f %f';
    case 'truncnorm', num = zeros(nt,6); strpat='%s %f %f %f %f %f %f';
    case 'yc1985',    num = zeros(nt,5); strpat='%s %f %f %f %f %f';
    case 'trexpub',   num = zeros(nt,7); strpat='%s %f %f %f %f %f %f %f';
end

for i=1:nt
    D=textscan(t{i},strpat);
    source{i}=D{1}{1};
    num(i,:)=horzcat(D{2:end});
end
mscl.source=source;
mscl.num=num;