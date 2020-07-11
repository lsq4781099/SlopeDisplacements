function[mulogIM,L]=run_DSHA(source,scen,h,opt)

r0        = gps2xyz(h.p,opt.ellipsoid);
Nscen     = size(scen,1);
[slist,~,sf]  = unique(scen(:,  1)');
[uscen,B] = unique(scen(:,1:2),'rows','stable');
Nsites    = size(r0,1);
IMs       = [opt.IM1;opt.IM2];
NIM       = length(IMs);
mulogIM   = zeros(Nscen,NIM,Nsites);
tau       = zeros(Nscen,NIM);
phi       = zeros(Nscen,NIM);


%source,M,normal,locxy,r0,VS30,opt

for i=slist
    ind = (sf==i);
    [mulogIM(ind,:,:),tau(ind,:),phi(ind,:)] = dsha_gmpe2(...
        source(i),scen(ind,2),scen(ind,3:5),scen(ind,6:8),r0,h,opt);
end

numgeom = vertcat(source.numgeom);
mec     = numgeom(sf,1); % list of mechanisms
L       = dsha_chol2(tau(B,:),phi(B,:),mec(B,:),uscen(:,2),r0,opt);

function [mulogIM,tau,phi] = dsha_gmpe2(source,M,rf,normal,r0,h,opt)

IMs           = [opt.IM1;opt.IM2];
NIM           = length(IMs);
Nsites        = size(r0,1);
Nscen         = size(M,1);
mulogIM       = zeros(Nscen,NIM,Nsites);
[tau,phi] = deal(nan(Nscen,NIM));

for i=1:Nsites
    source.media=h.value(i,:);
    if i==1, [mulogIM(:,:,i),tau,phi]=run_gmpe(source,M,normal,rf,r0(i,:),IMs,opt,h.param);
    else,    [mulogIM(:,:,i),~  ,~  ]=run_gmpe(source,M,normal,rf,r0(i,:),IMs,opt,h.param);
    end
end

function[mu,tau,phi]=run_gmpe(source,M,normal,rf,r0,IM,opt,hparam)
ellip = opt.ellipsoid;
Nscen = size(M,1);
NIM   = length(IM);
RA    = source.numgeom(3:4);
gmm   = source.gmm;
media = source.media;
%% ASSEMBLE GMPE PARAMERTER
switch source.obj
    case 1, param = param_circDSHA(r0,rf,M,normal,RA,gmm,media,ellip,hparam);
    case 2, param = param_circDSHA(r0,rf,M,normal,RA,gmm,media,ellip,hparam);
    case 3, param = param_circDSHA(r0,rf,M,normal,RA,gmm,media,ellip,hparam);
    case 4, param = param_circDSHA(r0,rf,M,normal,RA,gmm,media,ellip,hparam);
    case 5, param = param_rect(r0,source,ellip,hparam);
    case 6, param = param_circDSHA(r0,rf,M,normal,RA,gmm,media,ellip,hparam);
    case 7, param = param_circDSHA(r0,rf,M,normal,RA,gmm,media,ellip,hparam);
end

[mu,sig,tau,phi] = deal(zeros(Nscen,NIM));
for i=1:NIM
    [mu(:,i),sig(:,i),tau(:,i),phi(:,i)] = source.gmm.handle(IM(i),param{:});
end

std_exp   = 1;
sig_overw = 1;
sigma     = opt.Sigma;
if ~isempty(sigma)
    switch sigma{1}
        case 'overwrite', std_exp = 0; sig_overw = sigma{2};
    end
end
sig = sig.^std_exp*sig_overw;


%% This is an ASSUMPTION to deal with old GMMs that do not provide error separation in tau and phi
if any(isnan(tau(:))) || any(isnan(phi(:)))
    a   = 2/3;
    phi = sig/sqrt(1+a^2);
    tau = sqrt(sig.^2-phi.^2);
end

function L = dsha_chol2(TAU,PHI,MEC,MAG,r0,opt)

h           = computeh(r0);
funSpatial  = opt.Spatial;
funSpectral = opt.Spectral;

suTP = size(unique([TAU,PHI,MEC],'rows'),1);
if suTP==1
    switch func2str(funSpectral)
        case 'corr_none'             , Nunk = 1;
        case 'corr_BakerCornell2006' , Nunk = 1;
        case 'corr_BakerJayaram2008' , Nunk = 1;
        case 'corr_JayaramBaker2011' , Nunk = 1;
        case 'corr_Cimellaro2013'    , Nunk = 1;
        case 'corr_Abrahamson2014'   , Nunk = 1;
        case 'corr_BCHhydro2016'     , Nunk = 1;
        case 'corr_BakerBradley2017' , Nunk = 1;
        case 'corr_JaimesCandia2019' , Nunk = 1;
        case 'corr_Candia2019'       , Nunk = 1;
        case 'corr_GodaAtkinson2009' , Nunk = 1;
        case 'corr_Akkar2014'        , Nunk = 1;
        case 'corr_Ji2017'           , Nunk = length(MAG);
    end
else
    Nunk = length(MAG);
end

IMs    = [opt.IM1;opt.IM2];
[IMcorr1,IMcorr2] = meshgrid(IMs,IMs);
Nsites  = size(h,1);
NIM     = length(IMs);
L       = zeros(Nsites*NIM,Nsites*NIM,Nunk);

if all(isnan(TAU(:))) || all(isnan(PHI(:)))
    return
end
    
param(1:Nunk)= struct('opp',0,'residual','phi','direction','horizontal','isVS30clustered',1,'mechanism',[],'M',[]);
for i=1:Nunk
    switch MEC(i)
        case 1, param(i).mechanism = 'interface';
        case 2, param(i).mechanism = 'intraslab';
        case 3, param(i).mechanism = 'crustal';
    end
    param(i).M         = MAG(i);
end

%% Computes Spatial Correlation Structures
for jj=1:Nunk
    Lii    = zeros(Nsites,Nsites,NIM);
    tau    = TAU(jj,:);
    phi    = PHI(jj,:);
    paramj = param(jj);
    
    for i=1:NIM
        rho  = funSpatial(IMs(i), h, paramj);
        Cii  = (tau(i)^2+phi(i)^2*rho);
        if all(Cii==0)
            Lii(:,:,i)  = eye(size(Cii));
        else
            Lii(:,:,i)  = chol(Cii,'lower');
        end
    end
    
    % Computes Interperiod CorrelationStructures
    rhoIM=funSpectral(IMcorr1,IMcorr2,paramj);
    
    % Builds Extended Covariance Matrix
    C = zeros(Nsites*NIM,Nsites*NIM);
    for i=1:NIM
        II = (1:Nsites)+Nsites*(i-1);
        for j=1:NIM
            JJ       = (1:Nsites)+Nsites*(j-1);
            C(II,JJ) = rhoIM(i,j)*Lii(:,:,i)*Lii(:,:,j)';
        end
    end

    L(:,:,jj) = chol(C,'lower');
end

