function[area4]=process_obj6(area4,~)

%% process area4 sourcess
np    = size(area4.txt,1);
vert  = zeros(0,3);
vptr  = zeros(0,2);
mptr  = zeros(np,2);
xyzm  = zeros(0,3);
conn  = zeros(0,3);
aream = zeros(0,1);
hypm  = zeros(0,3);
center = zeros(np,3);

for j=1:np
    load(area4.adp{j},'geom')
    vertshift = size(vert,1);
    vert   = [vert;geom.vertices;geom.vertices(1,:);nan nan nan]; %#ok<*AGROW>
    center(j,:) =  mean(geom.vertices,1);
    vptr(j,:)=[1,size(geom.vertices,1)]+ vertshift;
    
    connshift = size(conn,1);
    conn      = [conn;geom.conn+size(xyzm,1)];
    xyzm      = [xyzm;geom.xyzm];
    aream     = [aream;geom.aream];
    hypm      = [hypm;geom.hypm];
    mptr(j,:) = [1,size(geom.conn,1)]+connshift;
end

area4.vptr  = vptr;
area4.mptr  = mptr;
area4.vert  = vert;
area4.xyzm  = xyzm;
area4.conn  = conn;
area4.aream = aream;
area4.hypm  = hypm;

% compute normal to each subarea (segments)
p1 = xyzm(conn(:,1),:);
p2 = xyzm(conn(:,2),:);
p3 = xyzm(conn(:,3),:);
n  = fastcross(p3-p1,p2-p1); % normal to each plane
nn = sum(n.^2,2).^0.5;
area4.normal = bsxfun(@rdivide,n,nn);
area4.center = center;
%area4.dip   = nan; %not yet implemented
