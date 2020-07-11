function[DSV]=dsv_line(xyzj,dip,ellipsoid)

xyz = 1/2*(xyzj(1,:)+xyzj(2,:));
x0 = xyz(1);
y0 = xyz(2);
z0 = xyz(3);

if ~isempty(ellipsoid.Code)
    a = ellipsoid.SemimajorAxis;
    b = ellipsoid.SemimajorAxis;
    c = ellipsoid.SemiminorAxis;
    
    E = [-y0/b^2,x0/a^2,0]'    ;E=E/norm(E);
    U = [x0/a^2,y0/b^2,z0/c^2]';U=U/norm(U);
    N = cross(U,E);
    ENU = [E,N,U];
else
    ENU = eye(3);
end

%% computes strike angle
n   = ENU(:,2);
rij = (xyzj(2,:)-xyzj(1,:))'/norm(xyzj(2,:)-xyzj(1,:));
strike =acosd(n'*rij); % because SeismicHazard uses ENU not END

k   = ENU(:,3);
K   = [0 -k(3) k(2);k(3) 0 -k(1);-k(2) k(1) 0];
ENU = (eye(3)+sind(-strike)*K+(1-cosd(-strike))*K^2)*ENU;

k = ENU(:,2);
K = [0 -k(3) k(2);k(3) 0 -k(1);-k(2) k(1) 0];
DSV = (eye(3)+sind(dip)*K+(1-cosd(dip))*K^2)*ENU;
