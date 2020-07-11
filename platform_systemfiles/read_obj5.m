function[area3]=read_obj5(t)

t=strrep(t,'area3 ','');
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
num   = zeros(Np,23);
vert  = cell(Np,1);
vptr  = zeros(Np,2);
node1 = 1;
for i=1:Np
     D=regexp(t{i},'\ ','split');
     txt{i}   = D{1};
     num(i,:) = str2double(D(2:end));
     
     lv      = reshape(num(i,12:23),3,4)';
     node2   = node1+4-1;
     vert{i} = [lv;lv(1,:);nan nan nan];
     vptr(i,:)=[node1 node2];
     node1   = node2+3;
end
area3      = createObj('source');
area3.txt  = txt;
area3.num  = num(:,1:11);
area3.vert = vertcat(vert{:});
area3.vptr = vptr;
