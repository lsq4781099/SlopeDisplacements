function[vout]=psda_J07Ia(ky,~,AIcm,varargin)
                          
% Jibson, R. W. (2007). Regression models for estimating coseismic 
% landslide displacement. Engineering geology, 91(2-4), 209-218.

%%
% Arias Intenity (AI) must be input in units of cm/s
% conversion of AI to m/s 
AI     = AIcm/100;
logEDP = 2.401*log10(AI)-3.481*log10(ky)-3.230;  % eq (9) Jibson (2007)
logsig = 0.656*ones(size(logEDP));

% log base convertion
lnEDP  = logEDP * log(10); 
sig    = logsig * log(10);

if nargin==3
    vout=[lnEDP,sig];
    return
end

y    = varargin{1};
dist = varargin{2};
if strcmp(dist,'pdf')
    vout = lognpdf(y,lnEDP,sig);
elseif strcmp(dist,'cdf')
    vout = logncdf(y,lnEDP,sig);
elseif strcmp(dist,'ccdf')
    vout = 1-logncdf(y,lnEDP,sig);
end
