function[x,dP]=trlognpdf_psda(param,nX)

x    = param(1);
covx = param(2);
sig  = x*covx;
pd   = trlogndist(x,covx);

if numel(param)==3
    Nsta = param(3);
    if covx>0 && Nsta>=2
        %% try 1
        % constraints: mean(x) = x'*dP = mu, std(x) = sig, sum(dP)=1
        % x-regularization
        x  = (1:Nsta)';
        b  = sig/std(x);
        a  = param(1)-b*mean(x);
        x  = a+b*x;
        
        % P-regularization
        P  = pdf(pd,x);
        on = ones(Nsta,1);
        f  = x;
        ab = [x'*P x'*f;on'*P on'*f]\[param(1);1];
        dP = ab(1)*P+ab(2)*f;
        
    else
        x   = param(1);
        dP  = 1;
    end
else % random numers
    if covx==0
        a  = param(1);
        x  = a*ones(1,nX);
    else
        x = random(pd,1,nX);
        % regularization
        b  = sig/std(x);
        a  = param(1)-b*mean(x);
        x  = a+b*x;
    end
    dP = ones(1,nX)/nX;
end

function[pd]=trlogndist(x,covx)

s  = x*covx;   % standard deviation
v  = s^2;      % variance

% https://www.mathworks.com/help/stats/lognormal-distribution.html
mu    = log(x^2/sqrt(v+x^2));
sigma = sqrt(log(v/x^2+1));
pd    = makedist('lognormal',mu,sigma);





