function[point1]=read_obj1(t)

t=strrep(t,'point1 ','');
t=strrep(t,' interface '       ,' 1 ');
t=strrep(t,' intraslab '       ,' 2 ');
t=strrep(t,' crustal '         ,' 3 ');
t=strrep(t,' strike-slip '     ,' 0 ');
t=strrep(t,' normal '          ,' -90 ');
t=strrep(t,' normal-oblique '  ,' -45 ');
t=strrep(t,' reverse '         ,' 90 ');
t=strrep(t,' reverse-oblique ' ,' 45 ');
t=strrep(t,' null '            ,' 1 ');
t=strrep(t,' wc1994 '          ,' 2 ');
t=strrep(t,' ellsworth '       ,' 3 ');
t=strrep(t,' hanksbakun2001 '  ,' 4 ');
t=strrep(t,' somerville1999 '  ,' 5 ');
t=strrep(t,' wc1994r '         ,' 6 ');
t=strrep(t,' wc1994s '         ,' 7 ');
t=strrep(t,' strasser2010 '    ,' 8 ');

Np=length(t);
txt  = cell(Np,1);
num  = zeros(Np,7);
vert = zeros(Np,3);
vptr = zeros(Np,2);

for i=1:Np
    D=textscan(t{i},'%s %f %f %f %f %f %f %f %f %f %f');
    txt{i}=D{1}{1};
    num(i,:)=horzcat(D{2:8});
    vert(i,:)=horzcat(D{9:11});
    vptr(i,:)=[i i];

end
point1=createObj('source');
point1.txt  = txt;
point1.num  = num;
point1.vert = vert;
point1.vptr = vptr;
