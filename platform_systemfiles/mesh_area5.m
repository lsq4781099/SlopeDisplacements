function[p,conn,vol,hyp,edge,normal]=mesh_area5(pv,s3,lmax,nref,th,slices,ellip)

%% GENERATES QUADRATIC MESH OF THE 3D SOURCE GIVEN BY p
if size(pv,1)==4
    [p,conn,area] = mesh_quad(pv,lmax,nref);
else
    m       = mean(pv);
    N       = size(pv, 1);
    pv_mov  = pv-m;
    covar   = pv_mov'*pv_mov/N;
    [~,~,V] = svd(covar);
    pv_rot  = pv_mov*V;
    turbomesh = max(abs(pv_rot(:,3)))>0.03*max(mean(abs(pv_rot(:,1:2))));
    
    if turbomesh  && s3>1e-4
        [p,conn,area]=mesh_tria_turbo(pv,lmax,nref,ellip);
    else
        [p,conn]=pmesh2(pv_rot(:,1:2),lmax,nref);
        % Computes areas
        p  = bsxfun(@plus,m,[p,p(:,1)*0]*V');
        V1 = p(conn(:,1),:);
        V2 = p(conn(:,2),:);
        V3 = p(conn(:,3),:);
        area = 1/2*fastcross(V2-V1,V3-V1);
        area= sum(area.^2,2).^0.5;
        
        % % Deletes ill triangles
        Area_T    = sum(area);
        I=(area<1e-6*Area_T);
        conn(I,:) = [];
        
        % remove repeated nodes within l_tol
        l_tol    = 1e-6*sqrt(Area_T);
        [p,conn] = fixmesh(p,conn,l_tol);
        V1       = p(conn(:,1),:);
        V2       = p(conn(:,2),:);
        V3       = p(conn(:,3),:);
        area     = 1/2*fastcross(V2-V1,V3-V1);
        area     = sum(area.^2,2).^0.5;
        % computes hypocenter
        hyp = zeros(size(conn,1),3);
        for i=1:3
            hyp=hyp+1/3*p(conn(:,i),:);
        end
    end
end
Area_T   = sum(area);
l_tol    = 1e-6*sqrt(Area_T);
[p,conn] = fixmesh(p,conn,l_tol);

%% computes normal
p1 = p(conn(:,1),:);
p2 = p(conn(:,2),:);
p3 = p(conn(:,3),:);
n  = fastcross(p3-p1,p2-p1); % normal to each plane
nn = sum(n.^2,2).^0.5;
normal = mean(bsxfun(@rdivide,n,nn),1);

%% MESH EXTRUSION
slices = slices+1;
th_vec = linspace(-th/2,th/2,slices);
pt = zeros(0,3);
for i=1:length(th_vec)
    newpt = p+normal*th_vec(i);
    pt    = [pt;newpt];
end

face1=pv([1:end,1],:)+normal*th/2;
face2=pv([1:end,1],:)-normal*th/2;
edge = [face1;nan nan nan;face2;nan nan nan];
for i=1:size(face1,1)
    edge=[edge;face1(i,:);face2(i,:);nan nan nan];
end

DT   = delaunayTriangulation(pt);
p    = DT.Points;
conn = DT.ConnectivityList;
a    = p(conn(:,1),:);
b    = p(conn(:,2),:);
c    = p(conn(:,3),:);
d    = p(conn(:,4),:);
vol  = 1/6*abs(sum((a-d).*fastcross(b-d,c-d),2));
hyp  = (a+b+c+d)/4;

function [p,t,area]=pmesh2(pv,hmax,nref)
pv = [pv;pv(1,:)];

p=[];
for i=1:size(pv,1)-1
    pp=pv(i:i+1,:);
    L=sqrt(sum(diff(pp,[],1).^2,2));
    if L>hmax
        n=ceil(L/hmax);
        pp=interp1([0,1],pp,(0:n)/n);
    end
    p=[p;pp(1:end-1,:)];  %#ok<*AGROW>
end

while 1
    t=delaunay(p);
    t=removeoutsidetris(p,t,pv);
    
    area=triarea(p,t);
    [maxarea,ix]=max(area);
    if maxarea<hmax^2/2, break; end
    pc=circumcenter(p(t(ix,:),:));
    p(end+1,:)=pc;
end

for iref=1:nref
    p = [p;edgemidpoints(p,t)];
    p = unique(p,'rows');
    t=delaunayn(p);
    t=removeoutsidetris(p,t,pv);
end

function pmid=edgemidpoints(p,t)

pmid=[(p(t(:,1),:)+p(t(:,2),:))/2;
    (p(t(:,2),:)+p(t(:,3),:))/2;
    (p(t(:,3),:)+p(t(:,1),:))/2];
pmid=unique(pmid,'rows');

function a=triarea(p,t)

d12=p(t(:,2),:)-p(t(:,1),:);
d13=p(t(:,3),:)-p(t(:,1),:);
a=abs(d12(:,1).*d13(:,2)-d12(:,2).*d13(:,1))/2;

function t=removeoutsidetris(p,t,pv)

pmid=(p(t(:,1),:)+p(t(:,2),:)+p(t(:,3),:))/3;
isinside=inpolygon(pmid(:,1),pmid(:,2),pv(:,1),pv(:,2));
t=t(isinside,:);

function pc=circumcenter(p)

dp1=p(2,:)-p(1,:);
dp2=p(3,:)-p(1,:);

mid1=(p(1,:)+p(2,:))/2;
mid2=(p(1,:)+p(3,:))/2;

s=[-dp1(2),dp2(2);dp1(1),-dp2(1)]\(-mid1+mid2)';
pc=mid1+s(1)*[-dp1(2),dp1(1)];
