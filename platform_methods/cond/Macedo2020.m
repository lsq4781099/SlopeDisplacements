function [lny,sigma,tau,phi] = Macedo2020(To,M,~,Vs30,func,varargin)

st = dbstack;
[isadmisible,units] = isIMadmisible(To,st(1).name,[nan nan],[nan nan],[nan nan],[nan nan]);
if isadmisible==0
    lny   = nan(size(M));
    sigma = nan(size(M));
    tau   = nan(size(M));
    phi   = nan(size(M));
    return
end

%% Model Coefficients
c1  = -0.79771;
c2  = 0.57546;
c3  = 0.48913;
c4  = -0.30894;
c5  = 0.11076;
tau = 0.18282;
phi = 0.21353;
rho = 0.52;

%% Median GMPE model for CAV (Shallow Crustal)
[lnPGA,sPGA]  = func(0,varargin{:}); % evaluates GMM for PGA
[lnSA1,sSA1]  = func(1,varargin{:}); % evaluates GMM for Sa(T=1)

lny   = c1 + c2*lnPGA + c3*M + c4*log(Vs30) + c5*lnSA1;
sigma = sqrt(0.0784+c2^2*sPGA.^2+c5^2*sSA1.^2+0.1276*rho*sPGA.*sSA1);

% unit conversion
lny = lny+log(units);
