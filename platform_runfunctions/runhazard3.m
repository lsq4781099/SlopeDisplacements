function[deagg]=runhazard3(im,IM,h,opt,source,Nsource,site_selection,Ebin)

ellipsoid = opt.ellipsoid;
xyz       = gps2xyz(h.p,ellipsoid);
Nsite     = size(xyz,1);
NIM       = length(IM);
Nim       = size(im,1);
deagg     = cell(Nsite,Nim,NIM,Nsource);

e1 = Ebin(1,2)-0.5;
e2 = Ebin(end,1)+0.5;
de = (Ebin(2,2)-Ebin(2,1))/3;
sigma     = opt.Sigma;
epsilon   = linspace(e1,e2,round((e2-e1)/de+1));
emin      = epsilon (1);
emax      = epsilon (end);
deps      = mean(diff(epsilon));

ind  = zeros(Nsite,length(source));
for i=site_selection
    ind(i,:)=selectsource(opt.MaxDistance,xyz(i,:),source);
end

for k=site_selection
    ind_k      = ind(k,:);
    sptr       = find(ind_k);
    xyzk       = xyz(k,:);
    valuek     = h.value(k,:);
    for i=sptr
        source(i).media = valuek;
        deagg(k,:,:,i)  = runsourceDeagg(source(i),xyzk,IM,im,ellipsoid,sigma,emin,emax,deps,h.param);
    end
end


return

function[deagg]=runsourceDeagg(source,r0,T,im,ellip,sigma,emin,emax,deps,hparam)

gmm = source.gmm;

%% MAGNITUDE RATE OF EARTHQUAKES
Nper   = length(T);
Nim    = size(im,1);
NMmin  = source.NMmin;

%% ASSEMBLE GMPE PARAMERSER
gmpefun  = gmm.handle;

switch source.obj
    case 1, [param,rate_MR] = param_circ(r0,source,ellip,hparam);  % point1
    case 2, [param,rate_MR] = param_circ(r0,source,ellip,hparam);  % line1
    case 3, [param,rate_MR] = param_circ(r0,source,ellip,hparam);  % area1
    case 4, [param,rate_MR] = param_circ(r0,source,ellip,hparam);  % area2
    case 5, [param,rate_MR] = param_rect(r0,source,ellip,hparam);  % area3
    case 6, [param,rate_MR] = param_circ(r0,source,ellip,hparam);  % area4
    case 7, [param,rate_MR] = param_circ(r0,source,ellip,hparam);  % volume1
end

%% HAZARD INTEGRAL
std_exp   = 1;
sig_overw = 1;
PHI       = 0;
if ~isempty(sigma)
    switch sigma{1}
        case 'overwrite', std_exp = 0; sig_overw = sigma{2};
        case 'truncate' , emin = min([emin,sigma{2},emin+sigma{2}]); emax = sigma{2}; PHI = 0.5*(1-erf(sigma{2}/sqrt(2)));
    end
end

deagg    = cell(Nim,Nper);
ep1      = emin:deps:(emax-deps);
ep2      = (emin+deps):deps:emax;
epsilon  = 1/2*(ep1+ep2);
Neps     = length(epsilon);

ep1(1)   = -inf;
ep2(end) = inf;

if sig_overw == 1
    rate_e   = normcdf(ep2)-normcdf(ep1);
end

if sig_overw ~= 1
    rate_e = zeros(1,Neps);
    rate_e(epsilon==sig_overw)=1;
end

Zhyp = nan(size(param{1}));
switch gmm.type
    case 'regular'
    NMR      = size(param{1},1);
    Mag      = repmat(param{1},Neps,1);%creating copies for the magnitude deagregation
    Rup      = repmat(param{2},Neps,1);%creating copies for the distance deagregation
    if source.gmm.Rmetric(8)
        Zhyp = repmat(param{3},Neps,1);%creating copies for the depth deagregation
    end
    case 'udm'
    NMR      = size(param{5},1);
    Mag      = repmat(param{5},Neps,1);%creating copies for the magnitude deagregation
    Rup      = repmat(param{6},Neps,1);%creating copies for the distance deagregation
    case 'cond'
    NMR      = size(param{5},1);
    Mag      = repmat(param{5},Neps,1);%creating copies for the magnitude deagregation
    Rup      = repmat(param{6},Neps,1);%creating copies for the distance deagregation        
    case 'frn'
    NMR      = size(param{7}{1},1);
    Mag      = repmat(param{7}{1},Neps,1);%creating copies for the magnitude deagregation
    Rup      = repmat(param{7}{2},Neps,1);%creating copies for the distance deagregation            
end

epsilon  = repmat(epsilon',NMR,1); epsilon  = epsilon(:);
rate_e   = repmat(rate_e',NMR,1);  rate_e   = rate_e(:);
rate_MR  = repmat(rate_MR,Neps,1);
rate     = rate_MR.*rate_e;

for j=1:Nper
    [mu,std] = gmpefun(T(j),param{:});
    std      = std.^std_exp*sig_overw;
    mu       = repmat(mu,Neps,1);
    std      = repmat(std,Neps,1);
    lnSa     = mu+std.*epsilon;
    imj      = im(:,j); 
    
    for i=1:Nim
        im_i = log(imj(i));
        ccdf = (lnSa>=im_i)*1/(1-PHI);
        deagg{i,j} = [Mag,Rup,Zhyp,epsilon,NMmin*ccdf.*rate];%that will provide the rates for all possible combinations of M, D, and epsilon
    end
end

return

