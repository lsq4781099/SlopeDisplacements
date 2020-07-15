function[haz]=haz_PSDA_cdmM(handles)
handles.site_selection = 1:size(handles.h.p,1);
opt = handles.opt;

% finds what integration methods are declared in the logic tree
TCDM      = handles.tableCDM.Data(:,2:5);
[imethod,usedM]=getmethod(handles,TCDM);
haz=struct('imstandard',[],'IMstandard',[],'lambda',[],'deagg',[],'imvector',[],'IMvector',[],'corrlist',[],'MRD',[]);

if any(ismember(imethod,6))
    handles.opt        = opt_update(handles,usedM,opt);
    branch             = regexp(handles.tableCDM.Data(:,2),'\,','split');
    branch             = unique(str2double(vertcat(branch{:})),'rows','stable');
    handles.sys.branch = branch;
    handles.sys.weight = ones(size(branch,1),5);
        
    haz.imstandard = handles.opt.im;
    haz.IMstandard = handles.opt.IM;
    [haz.lambda,haz.deagg]=runlogictree2(handles.sys,handles.opt,handles.h,handles.site_selection);
    haz.IMstandard = handles.opt.IM;
end


function[imethod,usedM]=getmethod(handles,TCDM)

hazard  = regexp(TCDM(:,1),'\,','split');
hazard  = str2double(vertcat(hazard{:}));

ME      = handles.ME;
geomptr = unique(hazard(:,1));
s       = unique(vertcat(handles.sys.mech{geomptr}));
[~,B]   = intersect([1;2;3],s);
func    = {ME.str}';
Smodels = TCDM(:,B+1);
Smodels = unique(Smodels(:));
usedM   = cell(size(Smodels));

for i=1:length(usedM)
    [~,B]=intersect({handles.sys.SMLIB.id},Smodels{i});
    usedM{i}=handles.sys.SMLIB(B).str;
end

imethod = zeros(1,length(usedM));
for i=1:length(imethod)
    [~,B]=intersect(func,usedM{i});
    imethod(i)=ME(B).integrator;
end
usedM  = unique(usedM);
imethod = unique(imethod);

function[opt]=opt_update(handles,usedM,opt)
[~,B]    = intersect(handles.h.param,'Ts');
T2       = handles.h.value(:,B);
methods  = pshatoolbox_methods(5);

func = {methods.str}';
[~,b]=intersect(func,usedM);
IMfactor = zeros(0,1);
for i=1:length(b)
    IMfactor = [IMfactor;methods(b(i)).Safactor(:)];
end
IMfactor = unique(IMfactor);
Tnat     = unique(T2);
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
im  = nan(size(im0));
Nim = size(im,1);
for i=1:NIM
    ind = find(IM(i)==opt.IM);
    if ~isempty(ind)
        im(:,i)=im0(:,ind);
    else
        if IM(i)>=0
            im(:,i)=logsp(0.001,3,Nim);
        elseif IM(i)==-1
            im(:,i)=logsp(0.001,60,Nim);
        elseif IM(i)==-5
            im(:,i)=logsp(0.001,100,Nim);
        end
    end
end

opt.IM=IM;
opt.im=im;



