function [varargout] = set_I17(param,Sa1,CAVdp,varargin)
CAVdp   = CAVdp/980.66;  % cm/s to g*sec
Q       = param.Q;
B       = param.B;
LBS     = param.LBS;
HL      = param.HL;

if LBS<=16
    lnS   = -8.35+4.59*log(Q)-0.42*(log(Q))^2+0.072*LBS+0.58*log(tanh(HL/6))-0.02*B+0.84*log(CAVdp)+0.41*log(Sa1);
else
    lnS   = -7.48+4.59*log(Q)-0.42*(log(Q))^2+0.014*LBS+0.58*log(tanh(HL/6))-0.02*B+0.84*log(CAVdp)+0.41*log(Sa1);
end
sigma = 0.50;

if nargin==3
    varargout{1}=lnS;
    varargout{2}=sigma;
    return
end

sett  = varargin{1};
dist  = varargin{2};
Nd    = length(sett);

switch dist
    case 'convolute'
        vout  = zeros(1,Nd);
        MRD   = varargin{3};
        for k = 1:Nd
            xhat    = (log(sett(k))-lnS)./sigma;
            ccdf    = 0.5*(1-erf(xhat/sqrt(2)));
            dlambda = ccdf.*MRD;
            vout(k) = sum(dlambda(:));
        end
        varargout{1} = vout;
        
    case 'pdf'
        varargout{1} = lognpdf(sett,lnS,sigma);
    case 'cdf'
        varargout{1} = logncdf(sett,lnS,sigma);
    case 'ccdf'
        varargout{1} = 1-logncdf(sett,lnS,sigma);
end