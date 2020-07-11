function[mpdf,mcdf,Mo]=magtable(~,param)

Mo         = 0;
Mmin       = param(1);
binwidth   = param(2);
occurrance = param(3:end);
Nm         = length(occurrance);
m          = Mmin:binwidth:(Mmin+binwidth*(Nm-1));
ccdf       = occurrance/occurrance(1);
mcdf       = 1-ccdf;
mpdf       = timeDer22(m,mcdf,1);

function[u]=timeDer22(time,u,n)

time=time(:)';
u=u(:)';

Np = length(time);
h  = mean(diff(time));
e = ones(Np,1)/2;
A = spdiags([-e,e*0,e],[-1 0 1],Np,Np);
A(1,1:2)=[-1,1];
A(Np,Np-1:Np)=[-1,1];

for jj=1:n
    u   = (A*u')'/h;
end
