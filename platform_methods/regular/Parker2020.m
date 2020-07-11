function [lny,sigma,tau,phi]=Parker2020(To,M,Rrup,Zhyp,Vs30,Z10,Z25,mechanism,region,slab,aBackarc,aNankai,Mb,Basin,AleInSigma,EpiInSigma,Nsample)

%M          : Moment Magnitud
%Rrup       : Rupture distance [km]
%Event type : 'interface' or 'intraslab'
%Ztor       : Depth of the top of rupture [km]
%Vs30       : Average shear wave velocity in the top 30 meters [m/s]

st          = dbstack;
isadmisible = isIMadmisible(To,st(1).name,[0 10],[nan nan],[nan nan],[nan nan]);
if isadmisible==0
    lny   = nan(size(M));
    sigma = nan(size(M));
    tau   = nan(size(M));
    phi   = nan(size(M));
    return
end

To      = max(To,0.001);  % PGA is associated to To=0.01;
period  = [-1 0.001 0.01,0.02,0.03,0.05,0.075,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.75,1,1.5,2,3,4,5,7.5,10];
T_lo    = max(period(period<=To));
T_hi    = min(period(period>=To));
index   = find(abs((period - T_lo)) < 1e-6);    % Identify the period

if T_lo==T_hi
    [lny,sigma,tau,phi]=gmpe(index,M,Rrup,Zhyp,Vs30,Z10,Z25,mechanism,region,slab,aBackarc,aNankai,Mb,Basin,AleInSigma,EpiInSigma,Nsample);
else
    [lny_lo,~,tau_lo,phi_lo] = gmpe(index,   M,Rrup,Zhyp,Vs30,Z10,Z25,mechanism,region,slab,aBackarc,aNankai,Mb,Basin,AleInSigma,EpiInSigma,Nsample);
    [lny_hi,~,tau_hi,phi_hi] = gmpe(index+1, M,Rrup,Zhyp,Vs30,Z10,Z25,mechanism,region,slab,aBackarc,aNankai,Mb,Basin,AleInSigma,EpiInSigma,Nsample);
    x          = log([T_lo;T_hi]);
    Y_sa       = [lny_lo,lny_hi]';
    Y_sigma2   = [tau_lo,tau_hi]';
    Y_sigma1   = [phi_lo,phi_hi]';
    lny        = interp1(x,Y_sa,log(To))';
    phi        = interp1(x,Y_sigma2,log(To))';
    tau        = interp1(x,Y_sigma1,log(To))';
    sigma      = sqrt(phi.^2+tau.^2);
end

function[lny,sigma,tau,phi]=gmpe(index,M,Rrup,Ztor,Vs30,Z10,Z25,mechanism,region,slab,aBackarc,aNankai,Mb,Basin,AleInSigma,EpiInSigma,Nsample)
[lny,sigma,tau,phi]=deal(M*0);
% 
% if strcmp(Mb,'default')
%     Mb = getMbDefault(mechanism, region, slab);
% end
% 
% switch Basin
%     case 'NoBasin'  ,Seattle_Basin = 0;
%     case 'InSeattle',Seattle_Basin = 1;
% end
% 
% lny = KBCG20_medPSA(index, M, Rrup, aBackarc, aNankai, Ztor, mechanism, Vs30, Z10, Z25, Mb, region, Seattle_Basin);
% 
% 
% function Mb = getMbDefault(mechanism, region, slab)
% 
% slab_saturation_regions_SBZ       = {'Aleutian', 'Alaska', '-999', 'Cascadia', 'Central_America_S', 'Central_America_N', 'Japan_Pac', 'Japan_Phi', 'New_Zealand_N', 'New_Zealand_S', 'South_America_N', 'South_America_S', 'Taiwan_W', 'Taiwan_E'};
% SaturationRegionList              = {'Alaska', 'Aleutian', 'Cascadia', 'Central_America_N', 'Central_America_S', 'Japan_Pac', 'Japan_Phi', 'Northern_Mariana', 'New_Zealand_N', 'New_Zealand_S', 'South_America_N', 'South_America_S', 'Taiwan_E', 'Taiwan_W'};
% SaturationRegionListDbRegionArray = {'1_Alaska', '1_Alaska', '2_Cascadia', '3_CentralAmerica&Mexico', '3_CentralAmerica&Mexico', '4_Japan', '4_Japan', '4_Japan', '5_NewZealand', '5_NewZealand', '6_SouthAmerica', '6_SouthAmerica', '7_Taiwan', '7_Taiwan'};
% 
% switch mechanism
%     case 'interface', Mbdef = 7.8; slab_saturation_regions_Mb = [8, 8, -999, 7.56, 7.5, 7.45, 8.31, 7.28, -999, -999, 8.45, 8.45, 8, 8];
%     case 'intraslab', Mbdef = 7.6; slab_saturation_regions_Mb = [7.98, 7.2, -999, 7.2, 7.6, 7.4, 7.65, 7.55, 7.6, 7.4, 7.3, 7.25, 7.7, 7.7];
% end

% switch region
%     case 'New-Zeeland',Mb=Mbdef;
%     case 'global'     ,Mb=Mbdef;
%     otherwise
%         [~,B] = intersect(SaturationRegionList,slab);
%         thisSaturationRegionDbRegion      = SaturationRegionListDbRegionArray{B};
%         if lower(thisSaturationRegionDbRegion) ~= lower(region)
%             Mb = Mbdef;
%         else
%             [~,B] = intersect(slab_saturation_regions_SBZ,slab);
%             Mb    = slab_saturation_regions_Mb(B);
%         end
%         if Mb==-999
%             Mb = Mbdef;
%         end
%         
% end

