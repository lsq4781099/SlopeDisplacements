function[strike,normal,dip]=dsv_area1(xyzj,ellipsoid)

% gpsj = xyz2gps(xyzj,ellipsoid);
p1   = xyzj(1,:);
p2   = xyzj(2,:);
p3   = xyzj(3,:);

normal = fastcross(p2-p1,p3-p2); % normal to each plane
normal = normal/norm(normal);
strike = (p2-p1)/norm(p2-p1);

if ~isempty(ellipsoid.Code)
    a = ellipsoid.SemimajorAxis;
    b = ellipsoid.SemimajorAxis;
    c = ellipsoid.SemiminorAxis;
    gpsj = xyz2gps(xyzj,ellipsoid)*diag([1 1 0]);
    xyzm = mean(gps2xyz(gpsj,ellipsoid),1);
    x0   = xyzm(1);
    y0   = xyzm(2);
    z0   = xyzm(3);
    E = [-y0/b^2,x0/a^2,0]'    ;E=E/norm(E);
    U = [x0/a^2,y0/b^2,z0/c^2]';U=U/norm(U);
    N = cross(U,E);
    ENU = [E,N,U];
else
    ENU = eye(3);
end
dip = acosd(normal*ENU(:,3));

if dip>90
    dip=180-dip;
end
