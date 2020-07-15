function[MRE,Cz]=runPCE(source,r0,IM,im,Nreal,ellip,hparam)

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
sqpi  = sqrt(pi); % squared root of pi
Nscen = length(param{1});
HSa   = [H(0,zrnd);H(1,zrnd);H(2,zrnd);H(3,zrnd);H(4,zrnd)];
Cz    = zeros(Nscen, Nim, 5);
MRE   = zeros(Nim,NIM,Nreal);

for j=1:NIM
    [lnY,sigma] = gmpefun(IM(j),param{:});
    mu   = mean(lnY,1)';
    smu  = std (lnY,0,1)';
    stot = sqrt(sigma.^2+smu.^2);
    a   = -smu.^2./(2*sigma.^2) - 1/2;
    for  i = 1:Nim
        lnz = log(im(i,j));
        b   = (lnz - mu) .* smu./(sigma.^2);
        c   = - (lnz - mu).^2 ./(2*sigma.^2);
        ee  = exp(c - b.^2./(4*a));
        Cz(:, i, 1) = 1/1  * (1-normcdf((lnz - mu)./stot));
        Cz(:, i, 2) = 1/1  * smu./(2*sigma*sqpi).* ee .*1./((-a).^0.5);
        Cz(:, i, 3) = 1/2  * smu./(2*sigma*sqpi).* ee .* b./(2*(-a).^1.5);
        Cz(:, i, 4) = 1/6  * smu./(2*sigma*sqpi).* ee .* ((-2*a.*(1 + 2*a) + b.^2)./(4*(-a).^2.5));
        Cz(:, i, 5) = 1/24 * smu./(2*sigma*sqpi).* ee .* (-b.*(6*a.*(1 + 2*a)-b.^2)./(8*(-a).^3.5));
    end
    
    Cz  = permute(nansum(bsxfun(@times,rate,Cz),1),[2 3 1]);
    mre = NMmin*Cz*HSa;
    ind = sum(mre<0,1)>0;
    mre(:,ind)=nan;
    MRE(:,j,:) = mre;
end

