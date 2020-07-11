function[p,conn,le,hyp]=mesh_line1(pv,lmax,nref)
% POLYLINE MESH OF A 3D LINE DEFINED BY 2 VERTICES PV(1,:) AND PV(2,:)

Nnodes = size(pv,1);
p    = zeros(0,3);

for i=1:Nnodes-1
    A=mesh_segment(pv(i:i+1,:),lmax,nref);
    if i==1
        p = [p;A]; %#ok<*AGROW>
    else
        p = [p;A(2:end,:)];
    end
end
conn   = [1:size(p,1)-1;2:size(p,1)]';
le     = sum((p(1:end-1,:)-p(2:end,:)).^2,2).^0.5;
hyp    = (p(1:end-1,:)+p(2:end,:))/2;


function[p,conn,le,hyp]=mesh_segment(pv,lmax,nref)
% LINEAR MESH OF A 3D LINE DEFINED BY 2 VERTICES PV(1,:) AND PV(2,:)

r1   = pv(1,:);
r2   = pv(2,:);
L    = norm(r2-r1);
LMAX = lmax/(nref+1);
nx   = ceil(L/LMAX)+1;
p    = zeros(nx,2);
p(:,1)=linspace(pv(1),pv(2),nx);
p(:,2)=linspace(pv(3),pv(4),nx);
p(:,3)=linspace(pv(5),pv(6),nx);

conn = [1:nx-1;2:nx]';

% computes Lengts
le = ones(nx,1)*L/(nx-1);

% computes hypocenter
hyp = zeros(nx-1,3);
for i=1:2
    hyp=hyp+1/2*p(conn(:,i),:);
end
