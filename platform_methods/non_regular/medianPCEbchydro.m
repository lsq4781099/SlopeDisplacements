function [lnY,sigma,tau,phi] = medianPCEbchydro(To,M,Rrup,VS30)

if  To<0 || To> 3
    lnY   = nan(size(M));
    sigma = nan(size(M));
    tau   = nan(size(M));
    phi   = nan(size(M));
    return
end

load PCE_bchydro_coefs.mat COEF
To         = max(To,0.01); %PGA is associated to To=0.01;
sigma      = 0.74*ones(size(M)); % this is period independent
period     = [0.01 0.02 0.05 0.075 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.6 0.75 1 1.5 2 2.5 3];
T_lo       = max(period(period<=To));
T_hi       = min(period(period>=To));
index      = find(abs((period - T_lo)) < 1e-6); % Identify the period
X1         = ones(size(M));
X2         = M;
X3         = (10-M).^2;
X4         = log(Rrup+10*exp(0.4*(M-6)));
X5         = M.*log(Rrup+10*exp(0.4*(M-6)));
X6         = Rrup;
X7         = log(VS30*ones(size(M))/760);
X          = [X1,X2,X3,X4,X5,X6,X7];

if T_lo==T_hi
    mean_theta = COEF(:,1,index);
    lnY        = X*mean_theta;
else
    x1    = log(T_lo);
    x2    = log(T_hi);
    x0    = log(To);
    S     = [1-(x0-x1)/(x2-x1),(x0-x1)/(x2-x1)];     
    mu1   = COEF(:,1,index);
    mu2   = COEF(:,1,index+1);
    mu    = [mu1,mu2]*S';
    lnY   = X*mu;
end

tau   = nan(size(M));
phi   = nan(size(M));
