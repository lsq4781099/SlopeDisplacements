function[lambda,deagg]=runlogictree2(sys,opt,h,sitelist)

%% variable initialization
IM        = opt.IM;
im        = opt.im;
Nsites    = size(h.p,1);
Nim       = size(im,1);
NIM       = length(IM);
Nbranch   = size(sys.branch,1);
weights   = sys.weight(:,5);
Nsource   = max(sum(sys.Nsrc,1));
if ~isfield(opt,'dflags')
    opt.dflags=[true true false];
end

%% do not run analysis if ind is empty
lambda  = nan (Nsites,Nim,NIM,Nsource,Nbranch);
deagg   = cell(Nsites,Nim,NIM,Nsource,Nbranch);

%% run logic tree
for i=1:Nbranch
    fprintf('%g\n',i/Nbranch)
    if weights(i)~=0
        source          = buildmodelin(sys,sys.branch(i,:),opt);
        deagg(:,:,:,:,i)= runhazard2(im,IM,h,opt,source,Nsource,sitelist);
    end
end

for i=1:numel(deagg)
    if ~isempty(deagg{i})
        lambda(i)=sum(deagg{i}(:,end));
    end
end

