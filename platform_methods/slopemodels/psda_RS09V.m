function[vout]=psda_RS09V(ky,~,pgv,pga,varargin)

% Rathje, E. M., & Saygili, G. (2009). Probabilistic assessment of 
% earthquake-induced sliding displacements of natural slopes. Bulletin 
% of the New Zealand Society for Earthquake Engineering, 42(1), 18.

% pgv must be input in cm/s
aratio = ky./pga;
%------------------------------------------------------------
e1 = -1.56;
e2 = -4.58;
e3 =-20.84;
e4 = 44.75;
e5 =-30.50;
e6 = -0.64;
e7 =  1.55;

lnd   = e1+e2*(aratio)+e3*(aratio).^2+e4*(aratio).^3+e5*(aratio).^4+e6*log(pga)+e7*log(pgv); 
std   = 0.405+0.524*aratio;
sig   = std.*ones(size(lnd));

if nargin==4
   vout=[lnd,sig];
   return
end

d    = varargin{1};
dist = varargin{2};
Nd   = length(d);

switch dist
    case 'convolute'
        vout  = zeros(1,Nd);
        MRD   = varargin{3};
        for k = 1:Nd
            xhat    = (log(d(k))-lnd)./sig;
            ccdf    = 0.5*(1-erf(xhat/sqrt(2)));
            dlambda = ccdf.*MRD.*(ky<=pga);
            vout(k)   = sum(dlambda(:));
        end
    case 'pdf'
        vout = lognpdf(d,lnd,sig);
    case 'cdf'
        vout = logncdf(d,lnd,sig);
    case 'ccdf'
        vout = 1-logncdf(d,lnd,sig);
end


