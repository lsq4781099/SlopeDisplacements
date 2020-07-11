function[area1]=process_obj3(area1,ellipsoid)

%% process area1 sourcess
np    = size(area1.txt,1);
mptr  = zeros(np,2);
xyzm  = zeros(0,3);
conn  = zeros(0,4);
aream = zeros(0,1);
hypm  = zeros(0,3);
center = zeros(np,3);

num   = area1.num;
vert  = area1.vert;
vptr  = area1.vptr;

for j=1:np
    lmax = num(j,6);
    nref = num(j,7);
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

area1.mptr  = mptr;
area1.xyzm  = xyzm;
area1.conn  = conn;
area1.aream = aream;
area1.hypm  = hypm;

% compute normal to each subarea (segments)
Nseg   = size(hypm,1);
strike = zeros(Nseg,3);
normal = zeros(Nseg,3);
dip    = zeros(Nseg,1);
for j=1:Nseg
    xyzj  = xyzm(conn(j,:),:);
    [strike(j,:),normal(j,:),dip(j)]=dsv_area1(xyzj,ellipsoid);
end
area1.strike = strike;
area1.normal = normal;
area1.dip    = dip;
area1.center = center;
