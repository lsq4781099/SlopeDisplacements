function[MRE]=runMCS(source,r0,IM,im,Nreal,ellip,hparam)

pd   = makedist('Normal');
t    = truncate(pd,-2,2);
zrnd = random(t,1,Nreal);

gmpe     = source.gmm;
NIM      = length(IM);
Nim      = size(im,1);
NMmin    = source.NMmin;
gmpefun  = gmpe.handle;

switch source.obj
    case 1, [param,rate] = param_circ(r0,source,ellip,hparam);  % point1
    case 2, [param,rate] = param_circ(r0,source,ellip,hparam);  % line1
    case 3, [param,rate] = param_circ(r0,source,ellip,hparam);  % area1
    case 4, [param,rate] = param_circ(r0,source,ellip,hparam);  % area2
    case 5, [param,rate] = param_rect(r0,source,ellip,hparam);  % area3
    case 6, [param,rate] = param_circ(r0,source,ellip,hparam);  % area4
    case 7, [param,rate] = param_circ(r0,source,ellip,hparam);  % volume1
end

%% HAZARD INTEGRAL
MRE   = zeros(Nim,NIM,Nreal);

for j=1:NIM
    [lnY,sigma] = gmpefun(IM(j),param{:});
    mu          = mean(lnY,1);
    smu         = std (lnY,0,1);
    lnzi        = log(im(:,j));
    
    MCS = zeros(Nim, Nreal);
    MU  = bsxfun(@plus,mu,zrnd'*smu)';
    for i = 1:Nim
        xhat = bsxfun(@times,lnzi(i)-MU,1./sigma);
        ccdf = 0.5*(1-erf(xhat/sqrt(2)));
        MCS(i,:) = NMmin* rate'*ccdf;
    end
    MRE(:,j,:)=MCS;
    
end

return







