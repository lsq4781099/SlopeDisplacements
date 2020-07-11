function[mpdf,mcdf,meanMo]=truncnorm(m,param)

Mmin    = param(1);
Mmax    = param(2);
Mchar   = param(3);
sigmaM  = param(4);

PHI  = normcdf([0,Mmin,Mmax],Mchar,sigmaM);
mpdf = normpdf(m,Mchar,sigmaM)/(PHI(3)-PHI(2));
mcdf = (normcdf(m,Mchar,sigmaM)-normcdf(Mmin,Mchar,sigmaM))/(PHI(3)-PHI(2));

% computation of moment rate meanMo  OJO aca
[mi,wi]= gaussquad(10,[0 Mmax]);
mpdfi  = normpdf(mi,Mchar,sigmaM)/(PHI(3)-PHI(2));
meanMo = wi'*(10.^(1.5*mi+16.05).*mpdfi);