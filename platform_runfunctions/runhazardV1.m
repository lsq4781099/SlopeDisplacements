function[lambda,MRE,MRD]=runhazardV1(im,IM,h,opt,source,Nsource,rho,varargin)

sigma      = opt.Sigma;
ellipsoid  = opt.ellipsoid;
xyz        = gps2xyz(h.p,ellipsoid);
Nim        = size(im,1);
rho        = min(rho,0.99999);
ind        = selectsource(opt.MaxDistance,xyz,source);
sptr       = find(ind);
lambda     = zeros(Nim,2,Nsource);
MRE        = nan(Nim,Nim,Nsource);
MRD        = nan(Nim,Nim,Nsource);

if nargin==7
    for i=sptr
        source(i).media = h.value;
        [lambda(:,:,i),MRE(:,:,i),MRD(:,:,i)]=runsource(source(i),xyz,IM,im,rho,sigma,ellipsoid,h.param);
    end
    lambda = nansum(lambda,3);
    MRE    = nansum(MRE   ,3);
    MRD    = nansum(MRD   ,3);
end

if nargin==8 % this option runs only the sources listed in varargin{1}, this is used in the PSDA GUI to pre-compute vector-valued hazard
    mechlist=varargin{1};
    for i=sptr
        if any(source(i).numgeom(1)==mechlist)
            source(i).media = h.value;
            [lambda(:,:,i),MRE(:,:,i),MRD(:,:,i)]=runsource(source(i),xyz,IM,im,rho,sigma,ellipsoid,h.param);
        end
    end
end

return

function[lambda,MRE,MRD]=runsource(source,r0,IM,im,rho,sigma,ellip,hparam)

gmpe        = source.gmm;
[Nim0,ncols]= size(im);
Nim         = 40;
NMmin       = source.NMmin;
gmpefun     = gmpe.handle;

%% ASSEMBLE GMPE PARAMERTER
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
[mu1,sig1] = gmpefun(IM(1),param{:});
[mu2,sig2] = gmpefun(IM(2),param{:});

if ncols==1
    im1 = im(:,1);
    im2 = im(:,1);
else
    im1 = im(:,1);
    im2 = im(:,2);
end

% this trick is to compute MRE accurately using the 2D Simpson's rule
% The goal is to evaluate MRE = int(int(MRD(x,y),x,im1,inf),y,im2,inf)
% ------------------------------------------------------------------------
% scalar seismic hazard curves
pd=makedist('normal',0,1);
nsig   = 6;
if ~isempty(sigma) && strcmp(sigma{1},'truncate')
    pd=makedist('normal',0,1);
    pd=truncate(pd,-inf,sigma{2});
    nsig   = sigma{2}+1;
end

if ~isempty(sigma) && strcmp(sigma{1},'overwrite')
    sig1=sig1*sigma{2};
    sig2=sig2*sigma{2};
end

maxIM1  = max(max(exp(mu1+nsig*sig1)),1e2*max(im1));
maxIM2  = max(max(exp(mu2+nsig*sig2)),1e2*max(im2));
im1     = logsp(im1(1),maxIM1,Nim)';
im2     = logsp(im2(1),maxIM2,Nim)';
Nim     = length(im1);
[X1,X2] = meshgrid(im1,im2);

mu1   = permute(mu1  ,[2 3 1]);
sig1  = permute(sig1 ,[2 3 1]);
mu2   = permute(mu2  ,[2 3 1]);
sig2  = permute(sig2 ,[2 3 1]);
rate  = permute(rate ,[2 3 1]);

if all(sig1==0) && all(sig2==0)
    on  = ones(Nim0,Nim0);
    MRD = zeros(Nim0,Nim0);
    MRE = zeros(Nim0,Nim0);
    Nscen = size(rate,3);
    for i=1:Nscen
        MREi = NMmin*rate(i)*on;
        MREi(:,log(im(:,1))>mu1(i))=0;
        MREi(log(im(:,2))>mu2(i),:)=0;
        MRE = MRE+MREi;
    end
else
    
    xhat  = (log(X1)-mu1)./sig1;
    fSa1  = 1./(X1.*sig1).*pdf(pd,xhat);
    mu21  = mu2 + rho*sig2.*xhat;          % eq 6
    sig21 = sig2.*sqrt(1-rho^2);           % eq 7
    xhat2 = (log(X2)-mu21)./sig21;         % eq 5
    fSa2  = 1./(X2.*sig21).*pdf(pd,xhat2);
    MRD   = NMmin*nansum(fSa1.*fSa2.*rate,3);  % eq 3
    MRE   = zeros(Nim,Nim);
    
    x1MRD = log(im1)';
    x2MRD = log(im2);
    F     = X1.*X2.*MRD;
    for i=1:Nim-1
        x2 = x2MRD(i:end);
        for j=1:Nim-1
            x1 = x1MRD(j:end);
            MRE(i,j) = trapz(x1,trapz(x2,F(i:end,j:end),1));
        end
    end
    [Xo1,Xo2] =meshgrid(im(:,1),im(:,2));
    MRD = interp2(log(X1),log(X2),MRD,log(Xo1),log(Xo2),'spline',nan);
    MRE = interp2(log(X1),log(X2),MRE,log(Xo1),log(Xo2),'spline',nan);
end

% scalar seismic hazard
lambda = zeros(Nim0,2);
for i=1:Nim0
    x       = im(i,1);
    xhat    = (log(x)-mu1)./sig1;
    ccdf1   = 1-cdf(pd,xhat);
    deagg1  = ccdf1.*rate;
    
    x       = im(i,2);
    xhat    = (log(x)-mu2)./sig2;
    ccdf2   = 1-cdf(pd,xhat);
    deagg2  = ccdf2.*rate;
    lambda(i,:) = NMmin*nansum([deagg1,deagg2],3);
end

return


