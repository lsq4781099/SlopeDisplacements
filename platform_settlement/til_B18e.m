function [LnTilt,sigma] = til_B18e(param,cav)

% input cav in cm/s

logSadj = set_B18(param,cav);
B       = param.B;
th1     = param.th1;

alpha1 =0.509;
alpha2 =-0.936;
alpha3 =-0.102;

% stdS    = 0.6746; % Settlement
% sigma_r	=0.287;
% rho_r_S	=0.194;
% DST=Ds(1)-Hs(1)/2;
% b26 is Ds
% b25 is Hs
%sigma = sqrt(sigma_r^2+alpha1^2*stdS^2+2*stdS*rho_r_S*sigma_r);

LnTilt = alpha1*logSadj+alpha2*log(B)+alpha3*th1;
sigma = 0.50;
