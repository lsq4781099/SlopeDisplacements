function runCMS_method3(handles)

% set up data
pkl         = handles.sys.weight(:,5);
pkd         = zeros(size(pkl));
ellip       = handles.opt.ellipsoid;
site_ptr    = handles.pop_site.Value;
site        = handles.h.p(site_ptr,:);

h           = handles.h;
h.id        = h.id(site_ptr);
h.p         = h.p(site_ptr,:);
h.value     = h.value(site_ptr,:);

r0          = gps2xyz(h.p,ellip);
Tcond       = str2double(handles.Cond_Period.String);
T           = UHSperiods(handles);
T           = sort(unique([T;Tcond]));
Tcond_ptr   = T==Tcond;
Tr          = str2double(handles.Ret_Period.String);
Nmodels     = size(handles.sys.branch,1);
Nper        = length(T);
muCMSk      = zeros(Nper,Nmodels);
SIGMAk      = zeros(Nper,Nmodels);

opt = handles.opt;
opt.IM = Tcond;
col    = find(handles.opt.IM>=0);
im1    = handles.opt.im(:,col(1));
opt.im  = im1;
MRE = runlogictree1(handles.sys,opt,h,1);
MRE = permute(MRE,[2 5 1 3 4]);

Sat = zeros(1,Nmodels);
for i=1:Nmodels
    Sat(i)=robustinterp(MRE(:,i),opt.im,1/Tr,'loglog');
end
SaT   = exp(log(Sat)*pkl);
col   = find(handles.opt.IM>=0);
im1   = repmat(handles.opt.im(:,col(1)),1,length(T));
opt2  = handles.opt;

for model_ptr  = 1:Nmodels
    branch  = handles.sys.branch(model_ptr,1:3);
    sources = buildmodelin(handles.sys,branch,handles.opt.ShearModulus);
    Nsource = length(sources);
    
    % compute Hazard Deagregation for T* at Return Period Tr
    im2     = Sat(model_ptr);
    opt2.im = im2;
    opt2.IM = Tcond;
    lambda2 = nan(1,1,1,Nsource,1);
    deagg2  = runhazard2(im2,Tcond,h,opt2,sources,Nsource,1);
    for i=1:numel(deagg2)
        if ~isempty(deagg2{i})
            lambda2(i)=sum(deagg2{i}(:,3));
        end
    end
    lambda2 = permute(lambda2,[2,1,3,4]);
    lambda2 = nansum(lambda2,4);
    [handles,MRscen] = run_func(handles,sources,lambda2,deagg2);
    
    M     = MRscen(1);
    pkd(model_ptr) = MRscen(3);
    ptr1  = MRscen(4);
    ptr2  = MRscen(5);
    
    % compute GMPE prediction
    source         = sources(ptr1);
    gmpefun        = source.gmm.handle;
    source.mscl    = [M 1];
    source.aream   = source.aream(ptr2,:);
    source.hypm    = source.hypm(ptr2,:);
    source.normal  = source.normal(ptr2,:);
    source.media   = h.value;
    switch source.obj
        case 1, param = param_circ(r0,source,ellip,h.param);  % point1
        case 2, param = param_circ(r0,source,ellip,h.param);  % line1
        case 3, param = param_circ(r0,source,ellip,h.param);  % area1
        case 4, param = param_circ(r0,source,ellip,h.param);  % area2
        case 5, param = param_rect(r0,source,ellip,h.param);  % area3
        case 6, param = param_circ(r0,source,ellip,h.param);  % area4
        case 7, param = param_circ(r0,source,ellip,h.param);  % volume1
    end
    
    mu  = zeros(Nper,1);
    sig = zeros(Nper,1);
    for j=1:Nper
        [mu(j),sig(j)] = gmpefun(T(j),param{:});
    end
    
    % Step 4: compute e*
    eTs = (log(SaT)-mu(Tcond_ptr))/sig(Tcond_ptr);
    
    % Step 5: compute ebar
    rho     = zeros(size(T));
    methods = pshatoolbox_methods(4,handles.spatial_model.Value);
    func    = methods.func;
    Cond_param.opp       = 0;
    Cond_param.mechanism = 'interface';
    Cond_param.M         = M;
    Cond_param.residual  = 'phi';
    Cond_param.direction = 'horizontal';
    for i=1:length(rho)
        rho(i)  = func(Tcond,T(i),Cond_param);
    end
    
    % Step 6: compute and stores cms
    muCMSk(:,model_ptr)  = mu+eTs*rho.*sig;
    SIGMAk(:,model_ptr)  = sig.*sqrt(1-rho.^2);
end

% COMPUTE MEAN CS USING DEAGGREGATION WEIGHTS
pkd     = pkd/sum(pkd);
lncms   = muCMSk*pkd;
sigCMS  = sqrt((SIGMAk.^2+(bsxfun(@minus,muCMSk,lncms).^2))*pkd);
epsCMS  = 1.0; % can be changed to something else;
lncmss1 = lncms+epsCMS*sigCMS;
lncmss2 = lncms-epsCMS*sigCMS;
cms     = exp(lncms);
cmsk    = exp(muCMSk);
cmss1   = exp(lncmss1);
cmss2   = exp(lncmss2);

%  graphics

% ax1
delete(findall(handles.ax1,'tag','haz'));
handles.ax1.ColorOrderIndex=1;
plot(handles.ax1,im1(:,Tcond_ptr),MRE,'.-','tag','haz');
plot(handles.ax1,handles.ax1.XLim,[1/Tr 1/Tr],'k--','tag','haz')
Nmodels = size(handles.sys.branch,1);
L=legend(handles.ax1,compose('Branch %i',1:Nmodels));
L.Location = 'SouthWest';
L.Box      = 'off';
xlabel(handles.ax1,'Sa(T*)','fontsize',8)
ylabel(handles.ax1,'\lambda Sa(T*)','fontsize',8)

% ax2
delete(findall(handles.ax2,'tag','cms'));
delete(findall(handles.ax2,'tag','cmss'));
plot(handles.ax2,T  ,cmsk,'-','color',[1 1 1]*0.76,'tag','cms','handlevisibility','off')
handles.ax2.ColorOrderIndex=1;
plot(handles.ax2,T  ,cms,'.-','tag','cms')
plot(handles.ax2,[T;nan;T],[cmss1;nan;cmss2],'.-','tag','cmss')
L=legend(handles.ax2,'CMS (\mu)','CMS(\mu\pm\sigma)');L.Box='off';

ch=findall(handles.ax4,'tag','rho');  ch.XData = T;         ch.YData = rho;

set(handles.ax3,'fontsize',8,'visible','off')
ch=get(handles.ax3,'children');
set(ch,'Visible','off');

% uicontext
cF=get(0,'format');
format long g
num = [im1(:,Tcond_ptr),MRE];
data = num2cell([zeros(1,1+Nmodels);num]);
data(1,:)=[sprintf('Sa(T*=%g)',Tcond),compose('Branch %i',1:Nmodels)];
c = uicontextmenu;
uimenu(c,'Label','Copy data','Callback',{@data2clipboard_uimenu,data});
uimenu(c,'Label','Undock','Callback'   ,{@figure2clipboard_uimenu,handles.ax1});
set(handles.ax1,'uicontextmenu',c);

num = [T,cms,cmss1,cmss2];
data = num2cell([zeros(1,4);num]);
data(1,:) = {'T(s)','CMS','CMS+SIG','CMS-SIG'};
c = uicontextmenu;
uimenu(c,'Label','Copy data','Callback',{@data2clipboard_uimenu,data});
uimenu(c,'Label','Undock','Callback'   ,{@figure2clipboard_uimenu,handles.ax2});
set(handles.ax2,'uicontextmenu',c);

num = [T,rho];
data = num2cell([zeros(1,2);num]);
data(1,:)={'T(s)',sprintf('rho(T,T*)')};
c = uicontextmenu;
uimenu(c,'Label','Copy data','Callback',{@data2clipboard_uimenu,data});
uimenu(c,'Label','Undock','Callback'   ,{@figure2clipboard_uimenu,handles.ax4});
set(handles.ax4,'uicontextmenu',c);

format(cF);
end

function [handles,MRscen,rmin,rmax,Rcenter,Mcenter]=run_func(handles,source,lambda2,deagg2)

deagg{1}    = vertcat(deagg2{1,1,1,:});
indsource   = zeros(0,2);
Nsources    = size(deagg2,4);
for i=1:Nsources
    dg         = deagg2{1,1,1,i};
    if ~isempty(dg)
        Nscen      = size(dg,1);
        NM         = length(source(i).mscl(:,1));
        NR         = size(source(i).aream,1);
        ptr2       = sort(repmat((1:NR)',NM,1));
        indsource  = [indsource;[i*ones(Nscen,1),ptr2]];  %#ok<*AGROW>
    end
end

if isempty(deagg)
    return
end

% build deaggregation chart 'dchart'
rmin      = handles.Rbin(1);
rmax      = handles.Rbin(end);
Rcenter   = mean(handles.Rbin,2);
Mcenter   = mean(handles.Mbin,2);
handles.dchart = deagghazard(deagg,lambda2,Mcenter,Rcenter);

Mbar = nansum(deagg{1}(:,1).*deagg{1}(:,3))/nansum(deagg{1}(:,3));
Rbar = nansum(deagg{1}(:,2).*deagg{1}(:,3))/nansum(deagg{1}(:,3));
ind1 = and(handles.Mbin(:,1)<Mbar,handles.Mbin(:,2)>Mbar);
ind2 = and(handles.Rbin(:,1)<Rbar,handles.Rbin(:,2)>Rbar);

MBIN = handles.Mbin(ind1,:);
RBIN = handles.Rbin(ind2,:);
ind  = (deagg{1}(:,1)>=MBIN(1)).*(deagg{1}(:,1)<=MBIN(2)).*(deagg{1}(:,2)>=RBIN(1)).*(deagg{1}(:,2)<=RBIN(2))==1;

DG     = [deagg{1}(ind,:),indsource(ind,:)];
mean1  = mean(DG(:,1:2),1);
std1   = std(DG(:,1:2),0,1);
disc1  = bsxfun(@times,bsxfun(@minus,DG(:,1:2)  ,mean1),1./std1);
disc2  = bsxfun(@times,bsxfun(@minus,[Mbar,Rbar],mean1),1./std1);
DISC   = sum(bsxfun(@minus,disc1,disc2).^2,2);
[~,indmin] = min(DISC);
MRscen = DG(indmin,:);

end

