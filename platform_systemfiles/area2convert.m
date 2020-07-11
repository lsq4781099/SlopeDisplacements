function[vertices]=area2convert(gps,dip,usp,lsd,ellipsoid)

gps1 = gps; gps1(:,3)=-usp;
gps2 = gps; gps2(:,3)=-lsd;

if dip<=0
    error('Dip <= 0 not valid for area2 objects')
    vertices=[]; %#ok<*UNRCH>
    return
end

if dip==90
    if size(gps,1)==2
        vertices=[gps1(1,:);gps2;gps1(end:-1:2,:)];
    else
        vertices=[gps1;flipud(gps2)];
    end
    return
end
H =-gps1(1,3)+gps2(1,3);
W = H/sind(dip);
%% ROTATE DIP
Np   = size(gps1,1);
xyz1 = gps2xyz(gps1,ellipsoid);
xyz2 = gps2xyz(gps2,ellipsoid);

vrot = zeros(Np,3);
for i=1:Np
    v =-xyz2(i,:)+xyz1(i,:);
    v = v/norm(v);
    
    if i==1
        k = xyz1(2,:)-xyz1(1,:);
        k = k/norm(k);
    elseif and(i>1,i<Np)
        k = xyz1(i+1,:)-xyz1(i-1,:);
        k = k/norm(k);
    else
        k = xyz1(Np,:)-xyz1(Np-1,:);
        k = k/norm(k);
    end    
   
    kxv  = cross(k,v);
    kdv  = dot(k,v);
    vrot(i,:) = (v*cosd(dip-90)+kxv*sind(dip-90)+k*kdv*(1-cosd(dip-90)));
    
end

vrot=mean(vrot,1);

GPSL=xyz2gps(xyz1+0.8*W*vrot,ellipsoid);
GPSH=xyz2gps(xyz1+1.2*W*vrot,ellipsoid);
GPS=gps2*0;
for i=1:Np
    GPS(i,1)=interp1([GPSL(i,3),GPSH(i,3)],[GPSL(i,1),GPSH(i,1)],-lsd);
    GPS(i,2)=interp1([GPSL(i,3),GPSH(i,3)],[GPSL(i,2),GPSH(i,2)],-lsd);
    GPS(i,3)=-lsd;
end

if size(gps,1)==2
    vertices=[gps1(1,:);GPS;gps1(end:-1:2,:)];
else
    vertices=[gps1;flipud(GPS)];
end
