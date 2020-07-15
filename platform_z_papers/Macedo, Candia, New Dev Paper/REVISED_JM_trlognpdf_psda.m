% function[x,pdf,dP,xo,difxo]=trlognpdf_psda(param)
% 
% if all(param(2:3)>0)
%     
%     % truncated lognormal pdf
%     %mu    = log(param(1));
%     mu    = param(1);
%     sig   = param(2);
%     Nsta  = param(3);
%     param = logninv([0.0005 0.9995],mu,sig);
%     
%     % log normal pdf
%     xo = logsp(param(1),param(2),Nsta+1)';
%     difxo=xo(2:end)-xo(1:end-1);
%     x  = 1/2*(xo(1:end-1)+xo(2:end));
%     pdf  = lognpdf(x,mu,sig);
%     
%     area = logncdf(xo(end),mu,sig)-logncdf(xo(1),mu,sig);
%     dP   = logncdf(xo(2:end),mu,sig)-logncdf(xo(1:end-1),mu,sig);
%     
%     % area correction
%     pdf = pdf / area;
%     dP  = dP/area;
%     
% else
%     x   = param(1);
%     pdf = 1;
%     dP  = 1;
% end
