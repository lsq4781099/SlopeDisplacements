function[mpdf,mcdf,meanMo]=delta(M,Mchar)

zer   = zeros(size(M));
mpdf  = zer;
mcdf  = zer;

mpdf(M==Mchar)=1;
mcdf(M>=Mchar)=1;
meanMo = 10.^(1.5*Mchar+16.05);
