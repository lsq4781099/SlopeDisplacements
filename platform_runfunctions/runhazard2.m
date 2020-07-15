function[deagg]=runhazard2(im,IM,h,opt,source,Nsource,site_selection,OV)

xyz       = gps2xyz(h.p,opt.ellipsoid);
Nsite     = size(xyz,1);
NIM       = length(IM);
Nim       = size(im,1);
ind       = zeros(Nsite,length(source));

if isfield(opt,'dflag')
    dflag = opt.dflag;
else
    dflag = [true true false];
end

for i=site_selection
    ind(i,:)=selectsource(opt.MaxDistance,xyz(i,:),source);
end

deagg  = cell(Nsite,Nim,NIM,Nsource);
for k=site_selection
    ind_k      = ind(k,:);
    sptr       = find(ind_k);
    xyzk       = xyz(k,:);
    valuek     = h.value(k,:);
    for i=sptr
        source(i).media=valuek;
        deagg(k,:,:,i)=runsourceDeagg(source(i),xyzk,IM,im,opt.ellipsoid,opt.Sigma,h.param,dflag);
    end
end

return

function[deagg]=runsourceDeagg(source,r0,IM,im,ellip,sigma,hparam,dflag)

%% MAGNITUDE RATE OF EARTHQUAKES
NIM   = length(IM);
Nim   = size(im,1);
NMmin = source.NMmin;

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
deagg     = cell(Nim,NIM);
std_exp   = 1;
sig_overw = 1;
PHI       = 0;
if ~isempty(sigma)
    switch sigma{1}
        case 'overwrite', std_exp = 0; sig_overw = sigma{2};
        case 'truncate' , PHI = 0.5*(1-erf(sigma{2}/sqrt(2)));
    end
end

switch source.gmm.type
    case 'regular'
        Mag  = param{1};
        Rrup = param{2};
    case 'udm'
        Mag  = param{5};
        Rrup = param{6};
    case 'cond'
        Mag  = param{5};
        Rrup = param{6};
    case 'pce'
        Mag  = param{1};
        Rrup = param{2};        
    case 'frn'
        Mag  = param{6}{1};
        Rrup = param{6}{2};             
end

switch source.gmm.type
    case 'pce'
        for j=1:NIM
            switch func2str(source.gmm.handle)
                case 'PCE_nga'     ,[mu,sig] = medianPCE_nga(IM(j),param{:});
                case 'PCE_bchydro' ,[mu,sig] = medianPCE_bchydro(IM(j),param{:});
            end
            sig = sig.^std_exp*sig_overw;
            imj = im(:,j);
            for i=1:Nim
                xhat        = (log(imj(i))-mu)./sig;
                ccdf        = (0.5*(1-erf(xhat/sqrt(2)))-PHI)*1/(1-PHI);
                ccdf        = ccdf.*(ccdf>0);
                dsum        = [Mag,Rrup,mu];
                deagg{i,j}  = [dsum(:,dflag),NMmin*ccdf.*rate];
            end
        end
        
    otherwise
        for j=1:NIM
            [mu,sig] = source.gmm.handle(IM(j),param{:});
            sig = sig.^std_exp*sig_overw;
            imj = im(:,j);
            for i=1:Nim
                xhat        = (log(imj(i))-mu)./sig;
                ccdf        = (0.5*(1-erf(xhat/sqrt(2)))-PHI)*1/(1-PHI);
                ccdf        = ccdf.*(ccdf>0);
                dsum        = [Mag,Rrup,mu];
                deagg{i,j}  = [dsum(:,dflag),NMmin*ccdf.*rate];
            end
        end        
end

return

