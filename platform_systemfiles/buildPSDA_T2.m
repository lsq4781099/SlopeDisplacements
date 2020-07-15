function[T2]=buildPSDA_T2(paramPSDA,kyval,Tsval)

nTs         = paramPSDA.Tssamples;
nky         = paramPSDA.kysamples;
[Ts,dPTs]   = trlognpdf_psda([Tsval nTs]);
[ky,dPky]   = trlognpdf_psda([kyval nky]);
[ind1,ind2] = meshgrid(1:length(Ts),1:length(ky));
ind1        = ind1(:);
ind2        = ind2(:);
Nparam      = length(ind1);
T2          = [compose('set%i',1:Nparam)',num2cell([Ts(ind1),ky(ind2),dPTs(ind1).*dPky(ind2)])];