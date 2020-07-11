function[area4]=read_obj6(t)

t=strrep(t,'area4 ','');
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
adp   = cell(Np,1);
num   = zeros(Np,6);
for i=1:Np
     D=regexp(t{i},'\ ','split');
     txt{i} = D{1};
     adp{i} = D{7};
     num(i,:) = str2double(D([2:6,8]));
end
area4 = createObj('source');
area4.txt = txt;
area4.num = num;
area4.adp = adp;
