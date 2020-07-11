function [y]=robustinterp(X,Y,x,style)

ind = or(isnan(X+Y),isinf(X+Y));
X(ind)=[];
Y(ind)=[];

[Xu,ind]= unique(X,'stable');
Yu      = Y(ind);

switch style
    case 'loglog'
        zer=find(or(Xu==0,Yu==0));
        Xu(zer)=[];
        Yu(zer)=[];
        y = exp(interp1(log(Xu),log(Yu),log(x),'pchip'));
    case 'linear'
        if numel(Xu)==1
            y = NaN;
        else
            y = interp1(Xu,Yu,x);
        end
end
