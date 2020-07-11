function hd=LIBS_Mw_independent(fun,s,param,im,lambda)

Ns   = length(s);
hd   = zeros(1,Ns);

for i=1:Ns
    si    = s(i);
    PD    = fun(param,im,si,'ccdf');
    hd(i) = -trapz(lambda,PD);
end




