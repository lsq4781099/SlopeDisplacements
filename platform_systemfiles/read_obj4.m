function[area2]=read_obj4(t,ellipsoid)

t=strrep(t,'area2 ','');
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
%
Np    = length(t);
txt   = cell(Np,1);
num   = zeros(Np,10);
vert  = cell(Np,1);
vptr  = zeros(Np,2);
node1 = 1;
for i=1:Np
     D=regexp(t{i},'\ ','split');
     txt{i}  = D{1};
     lv      = str2double(D(12:end));
     Nnodes  = numel(lv)/2;
     
     num(i,:)= str2double(D(2:11));
     DIP     = num(i,6);
     USD     = num(i,7);
     LSD     = num(i,8);
     lv      = reshape(lv,2,Nnodes)';
     lv0     = area2convert(lv,DIP,USD,LSD,ellipsoid);
     vert{i} = [lv0;lv0(1,:);nan nan nan];
     node2   = node1+2*Nnodes-1;
     vptr(i,:)=[node1 node2];
     node1   = node2+3;
end

area2 = createObj('source');
area2.txt  = txt;
area2.num  = num;
area2.vert = vertcat(vert{:});
area2.vptr = vptr;