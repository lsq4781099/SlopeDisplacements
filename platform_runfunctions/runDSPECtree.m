function[SPEC,MRZ]=runDSPECtree(sys,opt,optdspec,h,sitelist)

%% variable initialization
p         = h.p;
Nsites    = size(p,1);
T0        = optdspec.primaryIM;
periods   = optdspec.periods;
epsilon   = optdspec.epsilon;
NIM       = length(periods);
Nbranch   = size(sys.branch,1);
branch    = sys.branch;
weights   = sys.weight(:,5);

%% run SPEC tree
home
Nsource = max(sum(sys.Nsrc,1));
SPEC = nan (Nsites,NIM,Nsource,Nbranch);
MRZ  = nan (Nsites,  7,Nsource,Nbranch);


if license('test','Distrib_Computing_Toolbox')
    for i=sys.isREG
        if weights(i)~=0
            sources = buildmodelin(sys,branch(i,:),opt,optdspec);
            [SPEC(:,:,:,i),MRZ(:,:,:,i)] = runDSPEC1(T0,periods,h,opt,sources,Nsource,epsilon,sitelist);
        end
    end
    
else
    for i=sys.isREG
        if weights(i)~=0
            sources = buildmodelin(sys,branch(i,:),opt);
            [SPEC(:,:,:,i),MRZ(:,:,:,i)] = runDSPEC1(T0,periods,h,opt,sources,Nsource,epsilon,sitelist);
        end
    end
end

end

