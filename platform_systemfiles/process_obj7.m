function[volume1]=process_obj7(volume1,ellipsoid)

%% process volume1 sourcess
np    = size(volume1.txt,1);
mptr  = zeros(np,2);
xyzm  = zeros(0,3);
conn  = zeros(0,4);
aream = zeros(0,1);
hypm  = zeros(0,3);
center = zeros(np,3);

num   = volume1.num;
vert  = volume1.vert;
vptr  = volume1.vptr;
normal = zeros(np,3);

for j=1:np
    lmax   = num(j,6);
    nref   = num(j,7);
    thick  = num(j,8);
    slices = num(j,9);
    ind    = vptr(j,1):vptr(j,2);
    pv     = gps2xyz(vert(ind,:),ellipsoid);
    center(j,:)    =  mean(vert(ind,:),1);
    connshift = size(xyzm,1);
    ptrshift  = size(conn,1);
    s3   = std(vert(ind,3));
    [xyzmj,connj,areamj,hypmj,~,normal(j,:)]=mesh_area5(pv,s3,lmax,nref,thick,slices,ellipsoid);
    xyzm = [xyzm;xyzmj]; %#ok<*AGROW>
    conn = [conn;connj+connshift];
    aream = [aream;areamj];
    hypm = [hypm;hypmj];
    mptr(j,:)=[1,size(connj,1)]+ptrshift;
end

volume1.mptr   = mptr;
volume1.xyzm   = xyzm;
volume1.conn   = conn;
volume1.aream  = aream;
volume1.hypm   = hypm;
volume1.normal = normal;
volume1.center = center;

