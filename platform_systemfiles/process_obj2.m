function[line1]=process_obj2(line1,ellipsoid)

%% process line sourcess
np    = size(line1.txt,1);
mptr  = zeros(np,2);
xyzm  = zeros(0,3);
conn  = zeros(0,2);
aream = zeros(0,1);
hypm  = zeros(0,3);
dip   = zeros(0,1);
center = zeros(np,3);

num   = line1.num;
vert  = line1.vert;
vptr  = line1.vptr;


for j=1:np
    lmax = num(j,7);
    nref = num(j,8);
    ind  = vptr(j,1):vptr(j,2);
    pv   = gps2xyz(vert(ind,:),ellipsoid);
    center(j,:) =  mean(vert(ind,:),1);
    connshift = size(xyzm,1);
    ptrshift  = size(conn,1);
    [xyzmj,connj,lej,hypmj]=mesh_line1(pv,lmax,nref);
    xyzm = [xyzm;xyzmj]; %#ok<*AGROW>
    conn = [conn;connj+connshift];
    aream   = [aream;lej];
    hypm = [hypm;hypmj];
    Nsegj = size(connj,1);
    dip  = [dip;ones(Nsegj,1)*num(j,6)];
    mptr(j,:)=[1,size(connj,1)]+ptrshift;
end

line1.mptr = mptr;
line1.xyzm = xyzm;
line1.conn = conn;
line1.aream= aream;
line1.hypm = hypm;

% compute dsv matrix for each subsegment
Nseg   = size(hypm,1);
strike = zeros(Nseg,3);
normal = zeros(Nseg,3);
for j=1:Nseg
    xyzj = xyzm(conn(j,:),:);
    dsv  = dsv_line(xyzj,dip(j),ellipsoid);
    strike(j,:) = dsv(:,2)';  % strike 
    normal(j,:) = dsv(:,3)';  % normal to the rupture plane
end
line1.strike = strike;
line1.normal = normal;
line1.dip    = dip;
line1.center = center;

