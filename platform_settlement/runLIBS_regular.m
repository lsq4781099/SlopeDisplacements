function[handles]=runLIBS_regular(handles)
delete(findall(handles.ax1,'type','line'));
drawnow

optlib = handles.optlib;
sett   = optlib.sett;
Nsett  = length(sett);
ShearMod  = handles.opt.ShearModulus;
IJK    = handles.IJK;
T1     = handles.T1;
T3     = handles.T3;

Nbranches = size(IJK,1);
setLIB    = handles.setLIB;
h         = handles.h;
Nsites    = size(h.p,1);
Nsources  = max(sum(handles.sys.Nsrc,1));
id        = {setLIB.label};
lambdaD   = nan(Nsites,Nsett,Nsources,Nbranches);


fprintf('\n');
spat  = 'Site %-17g | Branch %-3g of %-49g Runtime:  %-4.3f s\n';
t0 = tic;
fprintf('                                BUILDING SETTLEMENT HAZARD \n');
fprintf('-----------------------------------------------------------------------------------------------------------\n');

hd0 = zeros(size(sett));

% Computes Q and Sa Values for each scenario
nQ     = optlib.nQ;
nPGA   = optlib.nPGA;
nHL    = 2;
[II,JJ,KK]=meshgrid(1:nQ,1:nPGA,1:nHL);
II   = II(:);
JJ   = JJ(:);
KK   = KK(:);

QLBSHL = zeros(Nbranches,3,Nsites);

for site_ptr=1:Nsites
    brptr =1:Nbranches;
    brptr((cell2mat(T1(:,2))==0))=[];
    param = handles.param(site_ptr);
    Q     = trlognpdf_psda([param.Q nQ]);
    
    h           = handles.h;
    h.p         = h.p(site_ptr,:);
    h.param     = h.param;
    h.value     = h.value(site_ptr,:);
    
    % computes Magnitude and PGA percentiles
    pga_ptr = find(handles.haz.IM==0);
    lambda  = handles.haz.lambda(site_ptr,:,pga_ptr,:,:);
    lambda  = nansum(lambda,4);
    lambda  = permute(lambda,[2 5 1 3 4]);
    
    if isfield(handles.haz,'imvector')
        [LBS,HL1,HL2]=deal(zeros(nPGA,1));
        Tret = optlib.RetPeriod;
        for ll=1:nPGA
            opt0    = handles.opt;
            opt0.im = handles.haz.im(:,pga_ptr);
            opt0.IM = robustinterp(lambda(:,pga_ptr),opt0.im,1/Tret,'loglog');
            opt0.dflag = [true false true];
            source  = buildmodelin(handles.sys,handles.sys.branch(ll,:),opt0);
            deagg   = runhazard2(opt0.IM,0,h,opt0,source,Nsources,1);
            deagg   = vertcat(deagg{:});
            M      = (deagg(:,1)'*deagg(:,3))/sum(deagg(:,3));
            PGA    = exp((deagg(:,2)'*deagg(:,3))/sum(deagg(:,3)));
            [LBS(ll),HL1(ll),HL2(ll)]=calc_LBS_FS (param.CPT,param.wt, param.Df,M, PGA);
        end
        ind2 = IJK(:,2);
        QLBSHL(:,:,site_ptr) = [Q(II(ind2)) LBS(JJ(ind2)),HL1(JJ(ind2)).*(KK(ind2)==1)+HL2(JJ(ind2)).*(KK(ind2)==2)];
    else
        ind2 = IJK(:,2);
        QLBSHL(:,1,site_ptr) = Q(II(ind2));
    end
    
    for branch_ptr=brptr
        ti=tic;
        indT1    = IJK(branch_ptr,1); % pointer to hazaes model
        indT3    = IJK(branch_ptr,3); % pointer to settlement model
        geomptr  = handles.sys.branch(indT1,1);
        
        % site/building parameters
        B        = zeros(3,1);
        [~,B(1)] = intersect(id,T3{indT3,2}); fun1 = setLIB(B(1)).func; % interface
        [~,B(2)] = intersect(id,T3{indT3,3}); fun2 = setLIB(B(2)).func; % slab
        [~,B(3)] = intersect(id,T3{indT3,4}); fun3 = setLIB(B(3)).func; % crustal
        
        % run sources
        if ~isempty(handles.haz.lambda)
            indlist = nansum(permute(nansum(handles.haz.lambda(site_ptr,:,:,:,indT1),3),[2 4 1 3]),1);
            indlist = find(indlist);
        else
            indlist = nansum(permute(nansum(handles.haz.MRD(site_ptr,:,:,:,indT1),3),[2 4 1 3]),1);
            indlist = find(indlist);
        end
        
        for source_ptr = indlist
            mechptr = handles.sys.mech{geomptr}(source_ptr);
            Bs      = B(mechptr);
            switch mechptr
                case 1 , fun = fun1; % 'interface'
                case 2 , fun = fun2; % 'intraslab'
                case 3 , fun = fun3; % 'crustal'
            end
            
            integrator = setLIB(Bs).integrator;
            IMsite     = str2IM(setLIB(Bs).IM);
            hd         = hd0;
            
            if integrator==1 % magnitude independent models
                [~,IM_ptr] = intersect(handles.haz.IM,IMsite);
                im         = handles.haz.im(:,IM_ptr);
                lambda     = handles.haz.lambda(site_ptr,:,IM_ptr,source_ptr,indT1);
                lambda     = permute(lambda,[2,1]);
                
                if max(lambda)>0
                    param.Q = QLBSHL(branch_ptr,1,site_ptr);
                    hd      = LIBS_Mw_independent(fun,sett,param,im,lambda);
                end
            end
            
            if integrator==2
                cavdp   = handles.haz.imvector(:,1);
                sa1     = handles.haz.imvector(:,2);
                MRD     = handles.haz.MRD(site_ptr,:,:,source_ptr,indT1);
                MRD     = permute(MRD,[2 3 1]);
                
                if max(MRD(:))>0
                    param.Q   = QLBSHL(branch_ptr,1,site_ptr);
                    param.LBS = QLBSHL(branch_ptr,2,site_ptr);
                    param.HL  = QLBSHL(branch_ptr,3,site_ptr);
                    hd      = fun(param,sa1,cavdp,sett,'convolute',MRD);
                end
            end
            lambdaD(site_ptr,:,source_ptr,branch_ptr) = hd;
        end
        fprintf(spat,site_ptr,branch_ptr,Nbranches,toc(ti))
    end
end
handles.QLBSHL = QLBSHL;
handles.lambdaD=lambdaD;
fprintf('-----------------------------------------------------------------------------------------------------------\n');
fprintf('%-88sTotal:     %-4.3f s\n','',toc(t0));
handles.butt2.Value=1;
