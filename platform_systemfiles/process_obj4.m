function[area2]=process_obj4(area2,ellipsoid)

%% process area2 sourcess
np    = size(area2.txt,1);
mptr  = zeros(np,2);
xyzm  = zeros(0,3);
conn  = zeros(0,4);
aream = zeros(0,1);
hypm  = zeros(0,3);
center = zeros(np,3);

num   = area2.num;
vert  = area2.vert;
vptr  = area2.vptr;

for j=1:np
    lmax = num(j,9);
    nref = num(j,10);
    ind  = vptr(j,1):vptr(j,2);
    pv   = gps2xyz(vert(ind,:),ellipsoid);
    center(j,:) =  mean(vert(ind,:),1);
    connshift = size(xyzm,1);
    ptrshift  = size(conn,1);
    s3 = std(vert(ind,3));
    [xyzmj,connj,Aj,hypmj]=mesh_area1(pv,s3,lmax,nref,ellipsoid);
    xyzm = [xyzm;xyzmj]; %#ok<*AGROW>
    conn = [conn;connj+connshift];
    aream = [aream;Aj];
    hypm = [hypm;hypmj];
    mptr(j,:)=[1,size(connj,1)]+ptrshift;
end

area2.mptr   = mptr;
area2.xyzm   = xyzm;
area2.conn   = conn;
area2.aream  = aream;
area2.hypm   = hypm;
area2.center = center;

% compute normal to each subarea (segments)
Nseg   = size(hypm,1);
normal = zeros(Nseg,3);
dipt   = zeros(Nseg,1);
for j=1:Nseg
    xyzj  = xyzm(conn(j,:),:);
    [~,normal(j,:),dipt(j)]=dsv_area1(xyzj,ellipsoid);
end

area2.normal = normal;
area2.dip    = dipt;

