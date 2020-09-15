function[MRE,MREPCE]=runlogictree1(sys,opt,h,sitelist)

%% variable initialization
p         = h.p;
Nsites    = size(p,1);
IM        = opt.IM;
im        = opt.im;
Nim       = size(im,1);
NIM       = length(IM);
Nbranch   = size(sys.branch,1);
branch    = sys.branch;
weights   = sys.weight(:,5);

%% run logic tree
home
switch opt.SourceDeagg
    case 'off'
        Nsource=1;
    case 'on'
        Nsource=max(sum(sys.Nsrc,1));
        if Nsource>500
            str = sprintf('Model with %g seismic sources. Do you wish to save the source deaggregation',Nsource);
            answer = questdlg(str, ...
                'Source deaggregation', ...
                'No, I want my mommy','Yes, bring it on','No, I want my mommy');
            if strcmp(answer,'No, I want my mommy')
                Nsource=1;
                opt.SourceDeagg='off';
            end
        end
end
MRE     = nan (Nsites,Nim,NIM,Nsource,Nbranch);
MREPCE  = cell(1,Nbranch);

if license('test','Distrib_Computing_Toolbox')
    for i=sys.isREG
        if weights(i)~=0
            %fprintf('%i %i %i %i\n',branch(i,:));
            sources = buildmodelin(sys,branch(i,:),opt);
            MRE(:,:,:,:,i) = runhazard1(im,IM,h,opt,sources,Nsource,sitelist);
        end
    end
    
    for i=sys.isPCE
        if weights(i)~=0
            sources = buildmodelin(sys,branch(i,:),opt);
            MREPCE{i}=runhazard1PCE(im,IM,h,opt,sources,Nsource,sitelist);
        end
    end
else
    for i=sys.isREG
        if weights(i)~=0
            sources = buildmodelin(sys,branch(i,:),opt);
            MRE(:,:,:,:,i) = runhazard1(im,IM,h,opt,sources,Nsource,sitelist);
        end
    end
    
    for i=sys.isPCE
        if weights(i)~=0
            sources = buildmodelin(sys,branch(i,:),opt);
            MREPCE{i}=runhazard1PCE(im,IM,h,opt,sources,Nsource,sitelist);
        end
    end
    
end

end

