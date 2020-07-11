function[haz]=haz_LIBS(handles)
handles.site_selection = 1:size(handles.h.p,1);
opt = handles.opt;

fprintf('\n');
t0 = tic;
fprintf('                               PROBABILISTIC SEISMIC HAZARD\n');
fprintf('-----------------------------------------------------------------------------------------------------------\n');

T3       = handles.T3(:,2:4);
[imethod,usedM]  = getmethod(handles,T3);
haz=struct('im',[],'IM',[],'lambda',[]);

if any(imethod==1)
    handles.opt = opt_update(handles,usedM,opt);
    haz.im      = handles.opt.im;
    haz.IM      = handles.opt.IM;
    haz.lambda  = runlogictree1(handles.sys,handles.opt,handles.h,handles.site_selection);
end

if any(imethod==2) % Macedo & Bray 
    setLIB=handles.setLIB;
    handles.opt = opt_update(handles,usedM,opt);
    rhoCAVSA1   = handles.optlib.rhoCAVSA1;
    mechlist    = [];
    [~,B1]=intersect({setLIB.label},handles.T3(:,2)); D1=intersect({setLIB(B1).str},{'set_I17'});
    [~,B2]=intersect({setLIB.label},handles.T3(:,3)); D2=intersect({setLIB(B2).str},{'set_I17'});
    [~,B3]=intersect({setLIB.label},handles.T3(:,4)); D3=intersect({setLIB(B3).str},{'set_I17'});
    if ~isempty(D1), mechlist=[mechlist,1];end % 'interface'
    if ~isempty(D2), mechlist=[mechlist,2];end % 'instraslab'
    if ~isempty(D3), mechlist=[mechlist,3];end % 'crustal'
    
    haz.imvector   = handles.opt.im(:,[1 3]);
    haz.IMvector   = handles.opt.IM([1,3]);
    handles.opt.im = haz.imvector;
    handles.opt.IM = haz.IMvector;
    haz.corrlist = rhoCAVSA1;
    haz.MRD =runlogictree2V(handles.sys,handles.opt,handles.h,rhoCAVSA1,mechlist);
end

fprintf('-----------------------------------------------------------------------------------------------------------\n');
fprintf('%-88sTotal:     %-4.3f s\n','',toc(t0));

function[imethod,usedM]=getmethod(handles,T3)

methods = handles.setLIB;
geomptr = unique(handles.sys.branch(handles.IJK(:,1),1));
s       = unique(vertcat(handles.sys.mech{geomptr}));
[~,B]   = intersect([1;2;3],s);
func    = {methods.str}';
Smodels = T3(:,B);
Smodels = unique(Smodels(:));
usedM   = cell(size(Smodels));

for i=1:length(usedM)
    [~,B]=intersect({methods.label},Smodels{i});
    usedM{i}=methods(B).str;
end

imethod = zeros(1,length(usedM));
for i=1:length(imethod)
    [~,B]=intersect(func,usedM{i});
    imethod(i)=methods(B).integrator;
end

imethod = unique(imethod);

function[opt]=opt_update(handles,usedM,opt)

methods  = handles.setLIB;
func     = {methods.str}';
[~,b]=intersect(func,usedM);
IM = zeros(0,1);
for i=1:length(b)
    IMs = methods(b(i)).IM;
    IMs = regexp(IMs,'\-','split')';
    IM  = [IM;str2IM(IMs)]; %#ok<AGROW>
end
IM  = unique(IM);
NIM = length(IM);
im0 = opt.im;
Nim = size(im0,1);

im  = nan(Nim,NIM);

for i=1:NIM
    ind = find(IM(i)==opt.IM);
    if ~isempty(ind)
        im(:,i)=im0(:,ind);
    else
        if     IM(i)>=0  , im(:,i)=logsp(0.001,3,Nim);   % SA  Values
        elseif IM(i)==-4 , im(:,i)=logsp(1,2000,Nim);    % CAV Values
        end
    end
end

opt.IM=IM;
opt.im=im;



