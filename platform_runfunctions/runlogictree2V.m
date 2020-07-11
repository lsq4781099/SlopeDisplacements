function[MRD]=runlogictree2V(sys,opt,h,rho,mechlist)

Nsites   = size(h.p,1);
Nim      = size(opt.im,1);
Nbranch  = size(sys.branch,1);
Nsource  = max(sum(sys.Nsrc,1));
MRD      = zeros(Nsites,Nim,Nim,Nsource,Nbranch);


for k=1:Nbranch
    fprintf('%g\n',k/Nbranch)
    source  = buildmodelin(sys,sys.branch(k,:),opt);
    for i=1:Nsites
        hi.p     = h.p(i,:);
        hi.param = h.param;
        hi.value = h.value(i,:);
        [~,~,MRD(i,:,:,:,k)]=runhazardV1(opt.im,opt.IM,hi,opt,source,Nsource,rho,mechlist);
    end
end
