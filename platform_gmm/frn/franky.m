function [mu,sig]=franky(T,varargin)

Ndepend  = length(varargin)/4;
ind3     = (1:Ndepend)+2*Ndepend;
iptrs    = getIMptrs(T,varargin(ind3));
gmpefun  = varargin{iptrs};
ind2     = iptrs+1*Ndepend;
ind4     = iptrs+3*Ndepend;
param    = varargin{ind4};
[mu,sig] = gmpefun(T,param{:});


% unit conversion done indide franky
gmmunits = varargin{ind2};
mu       = mu-log(gmmunits);

