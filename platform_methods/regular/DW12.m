function[lny,sigma,tau,phi]=DW12(To,M,Rrup,media,SOF)

% Du, W. and Wang, G. (2013), A simple ground?motion prediction model 
% for cumulative absolute velocity and model validation. 
% Earthquake Engng Struct. Dyn., 42: 1189-1202. 
% DOI: https://doi.org/10.1002/eqe.2266

st = dbstack;
[isadmisible,units] = isIMadmisible(To,st(1).name,[nan nan],[nan nan],[nan nan],[nan nan]);
if isadmisible==0
    lny   = nan(size(M));
    sigma = nan(size(M));
    tau   = nan(size(M));
    phi   = nan(size(M));
    return
end

% model coefficients
c1 =  1.826;
c2 = -0.130;
c3 = -1.403;
c4 =  0.098;
c5 =  0.286;
h  =  8.455;
c6 =  0.481;
c7 = -0.155;
c8 =  0.095;


switch SOF
    case 'strike-slip',     FR = 0;   FN=0;
    case 'normal',          FR = 0;   FN=1;
    case 'normal-oblique',  FR = 0;   FN=1;
    case 'reverse',         FR = 1;   FN=0;
    case 'reverse-oblique', FR = 0.5; FN=0;
    case 'unspecified',     FR = 0;   FN = 0; % assumption
end


switch lower(media) % SGS site category 
    case {1,'b'}, SC = 0; SD = 0; 
    case {2,'c'}, SC = 1; SD = 0; 
    case {3,'d'}, SC = 0; SD = 1;  
end

% median model in units of gs
lny = c1 + c2*(8.5-M).^2 + (c3+c4*M).*log(sqrt(Rrup.^2+h^2)) + c5*SC + c6*SD + c7*FN + c8*FR;

% standard deviation model
CAV   = exp(lny);
sigma = 0.45-0.042*log(CAV/0.02);
sigma(CAV<=0.15) = 0.45;
sigma(CAV>=1.00) = 0.37;

tau = 0.247;
phi = sqrt(sigma.^2-tau^2);

% unit convertion
lny  = lny+log(units);



