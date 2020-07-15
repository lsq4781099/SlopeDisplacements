function[haz]=haz_PSDA(handles)
handles.site_selection = 1:size(handles.h.p,1);
opt = handles.opt;

fprintf('\n');
t0 = tic;
fprintf('                               PROBABILISTIC SEISMIC HAZARD\n');
fprintf('-----------------------------------------------------------------------------------------------------------\n');

% finds what integration methods are declared in the logic tree
T3       = handles.T3(:,2:4);
[imethod,usedM]  = getmethod(handles,T3);
haz=struct('imstandard',[],'IMstandard',[],'lambda',[],'deagg',[],'imvector',[],'IMvector',[],'corrlist',[],'MRD',[]);

if any(ismember(imethod,[1 2]))
    handles.opt    = opt_update(handles,usedM,opt,1);
    haz.imstandard = handles.opt.im;
    haz.IMstandard = handles.opt.IM;
    [haz.lambda,haz.deagg]=runlogictree2(handles.sys,handles.opt,handles.h,handles.site_selection);
    
end

if any(ismember(imethod,[3 4])) % Ellen´s rigid and flexible slopes
    SMLIB=handles.sys.SMLIB;
    handles.opt = opt_update(handles,usedM,opt,3);

    corrlist =[];
    T3models = unique(T3(:));
    for i=1:length(T3models)
        [~,B]=intersect({SMLIB.id},T3models{i});
        switch SMLIB(B).str
            case 'psda_RA2011F'
                corrlist =[corrlist;SMLIB(B).param.rho]; %#ok<*AGROW>
            case {'psda_RA2011R','psda_RS09V'}
                corrlist =[corrlist;SMLIB(B).param.rho];
        end
    end
    corrlist   = unique(corrlist);
    mechlist = [];
    [~,B1]=intersect({SMLIB.id},handles.T3(:,2)); D1=intersect({SMLIB(B1).str},{'psda_RA2011F','psda_RA2011R','psda_RS09V'});
    [~,B2]=intersect({SMLIB.id},handles.T3(:,3)); D2=intersect({SMLIB(B2).str},{'psda_RA2011F','psda_RA2011R','psda_RS09V'});
    [~,B3]=intersect({SMLIB.id},handles.T3(:,4)); D3=intersect({SMLIB(B3).str},{'psda_RA2011F','psda_RA2011R','psda_RS09V'});
    if ~isempty(D1), mechlist=[mechlist,1];end % 'interface'
    if ~isempty(D2), mechlist=[mechlist,2];end % 'instraslab'
    if ~isempty(D3), mechlist=[mechlist,3];end % 'crustal'
    
    haz.imvector = handles.opt.im;
    haz.IMvector = handles.opt.IM;
    haz.corrlist = corrlist;

    Nsites  = size(handles.h.p,1);
    Nim     = size(opt.im,1);
    Nbranch = size(handles.sys.branch,1);
    Nsource = max(sum(handles.sys.Nsrc,1));
    Ncorr   = length(corrlist);
    haz.MRD = zeros(Nsites,Nim,Nim,Nsource,Nbranch,Ncorr);
    for ii=1:Ncorr
        haz.MRD(:,:,:,:,:,ii)=runlogictree2V(handles.sys,handles.opt,handles.h,corrlist(ii),mechlist);
    end
end

fprintf('-----------------------------------------------------------------------------------------------------------\n');
fprintf('%-88sTotal:     %-4.3f s\n','',toc(t0));

function[imethod,usedM]=getmethod(handles,T3)

methods = pshatoolbox_methods(5);
geomptr = unique(handles.sys.branch(handles.IJK(:,1),1));
s       = unique(vertcat(handles.sys.mech{geomptr}));
[~,B]   = intersect([1;2;3],s);
func    = {methods.str}';
Smodels = T3(:,B);
Smodels = unique(Smodels(:));
usedM   = cell(size(Smodels));

for i=1:length(usedM)
    [~,B]=intersect({handles.sys.SMLIB.id},Smodels{i});
    usedM{i}=handles.sys.SMLIB(B).str;
end

imethod = zeros(1,length(usedM));
for i=1:length(imethod)
    [~,B]=intersect(func,usedM{i});
    imethod(i)=methods(B).integrator;
end

imethod = unique(imethod);

function[opt]=opt_update(handles,usedM,opt,mtype)
allTs  = handles.allTs;
nSites = size(allTs,1);
nTs    = handles.paramPSDA.Tssamples;
T2     = zeros(max(nTs,1),nSites);
for i=1:nSites
    T2(:,i)=trlognpdf_psda([allTs(i,:) nTs]);
end
Tnat = unique(T2(:));
methods  = pshatoolbox_methods(5);


switch mtype
    case 1, B = ismember(usedM,{'psda_BMT2017M','psda_BT2007','psda_BT2007M','psda_BM2019M','psda_J07M','psda_J07Ia','psda_RS09M','psda_AM1988'}); usedM = usedM(B);
    case 3, B = ismember(usedM,{'psda_RA2011F','psda_RA2011R','psda_RS09V'}); usedM = usedM(B);
end

func = {methods.str}';
[~,b]=intersect(func,usedM);
IMfactor = zeros(0,1);
for i=1:length(b)
    IMfactor = [IMfactor;methods(b(i)).Safactor(:)];
end
IMfactor = unique(IMfactor);
IM  = [];
for i=1:length(IMfactor)
    if IMfactor(i)<=0
        IM = [IM;IMfactor(i)]; %#ok<*AGROW>
    else
        IM = [IM;IMfactor(i)*Tnat]; %#ok<*AGROW>
    end
end

NIM = length(IM);
im0 = opt.im;
Nim = size(im0,1);

im  = nan(Nim,NIM);

for i=1:NIM
    ind = find(IM(i)==opt.IM);
    if ~isempty(ind)
        im(:,i)=im0(:,ind);
    else
        if IM(i)>=0
            im(:,i)=logsp(0.00001,3,Nim);
        elseif IM(i)==-1
            im(:,i)=logsp(0.001,2,Nim);
        elseif IM(i)==-5
            im(:,i)=logsp(0.001,100,Nim);
        end
    end
end

opt.IM=IM;
opt.im=im;



