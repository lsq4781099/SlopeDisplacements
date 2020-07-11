function runCMS_method1(handles)

% set up data
ellip       = handles.opt.ellipsoid;
site_ptr    = handles.pop_site.Value;

h           = handles.h;
h.id        = h.id(site_ptr);
h.p         = h.p(site_ptr,:);
h.value     = h.value(site_ptr,:);

r0          = gps2xyz(h.p,ellip);
model_ptr   = handles.pop_branch.Value;
Tcond       = str2double(handles.Cond_Period.String);
T           = UHSperiods(handles);
T           = sort(unique([T;Tcond]));
Tcond_ptr   = T==Tcond;
Tr          = str2double(handles.Ret_Period.String);
col         = find(handles.opt.IM>=0);
im1         = repmat(handles.opt.im(:,col(1)),1,length(T));
opt         = handles.opt;
branch      = handles.sys.branch(model_ptr,1:3);
geom_ptr    = branch(1);
Nsource     = sum(handles.sys.Nsrc(:,geom_ptr));

% compute seismic hazard for all Sa(T*)
sources = buildmodelin(handles.sys,branch,handles.opt.ShearModulus);
lambda1 = runhazard1(im1,T,h,opt,sources,Nsource,1);
lambda1 = permute(lambda1,[2 3 1]);

% compute Hazard Deagregation for T* at Return Period Tr
im2     = robustinterp(lambda1(:,Tcond_ptr),im1(:,Tcond_ptr),1/Tr,'loglog');
opt2    = opt;
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
[handles,MRscen,rmin,rmax,Rcenter,Mcenter] = run_func(handles,sources,lambda2,deagg2);

M     = MRscen(1);
rrup  = MRscen(2);
ptr1  = MRscen(4);
ptr2  = MRscen(5);

% compute UHS
uhs  = uhspectrum(handles.opt.im,lambda1,1/Tr);

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

mu  = zeros(size(uhs));
sig = zeros(size(uhs));
for j=1:length(T)
    [mu(j),sig(j)] = gmpefun(T(j),param{:});
end
gmpe    = exp(mu);

% Step 4: compute e*
eTs = (log(uhs(Tcond_ptr))-mu(Tcond_ptr))/sig(Tcond_ptr);

% Step 5: compute ebar
rho     = zeros(size(T));
methods = pshatoolbox_methods(4,handles.spatial_model.Value);
func    = methods.func;
Cond_param.opp       = 0;
Cond_param.mechanism = 'interface';%model.source(ptr1).mechanism;
Cond_param.M         = M;
Cond_param.residual  = 'phi';
Cond_param.direction = 'horizontal';
for i=1:length(rho)
    rho(i)  = func(Tcond,T(i),Cond_param);
end

% Step 6: compute cms
sigCMS  = sig.*sqrt(1-rho.^2);
epsCMS  = 1.0; % can be changed to something else;
lncms   = mu+eTs*rho.*sig;
lncmss1 = lncms+epsCMS*sigCMS;
lncmss2 = lncms-epsCMS*sigCMS;
cms     = exp(lncms);
cmss1   = exp(lncmss1);
cmss2   = exp(lncmss2);

%[T,cms,cmss1,cmss2];

%%  graphics

handles.control_R.String = sprintf('%4.3g',rrup);
handles.control_M.String = sprintf('%4.3g',M);

% ax1
delete(findall(handles.ax1,'tag','haz'));
handles.ax1.ColorOrderIndex=1;
plot(handles.ax1,im1(:,Tcond_ptr),lambda1(:,Tcond_ptr),'.-','tag','haz');

plot(handles.ax1,handles.ax1.XLim,[1/Tr 1/Tr],'k--','tag','haz')

xlabel(handles.ax1,'Sa(T*)','fontsize',8)
ylabel(handles.ax1,'\lambda Sa(T*)','fontsize',8)


% ax2
handles.ax2.ColorOrderIndex=1;
delete(findall(handles.ax2,'tag','cms'));
delete(findall(handles.ax2,'tag','cmss'));

plot(handles.ax2,T        ,cms,'.-','tag','cms')
plot(handles.ax2,[T;nan;T],[cmss1;nan;cmss2],'.-','tag','cmss')
L=legend(handles.ax2,'CMS (\mu)','CMS(\mu\pm\sigma)');L.Box='off';

ch=findall(handles.ax4,'tag','rho');  ch.XData = T;         ch.YData = rho;

set(handles.ax3,'fontsize',8,'visible','on')
handles.b=bar3(handles.ax3,Rcenter,handles.dchart(:,:,1));
set(handles.ax3,...
    'xtick',1:2:length(Mcenter),...
    'xticklabel',Mcenter(1:2:end),...
    'ytick',Rcenter(1:2:end),...
    'projection','perspective',...
    'fontsize',8,...
    'ylim',[rmin,rmax],'visible','on')
xlabel(handles.ax3,'Magnitude','fontsize',8);
ylabel(handles.ax3,'R (km)','fontsize',8);
zlabel(handles.ax3,'Hazzard Deagg','fontsize',8)

% uicontext
cF=get(0,'format');
format long g
num = [im1(:,Tcond_ptr),lambda1(:,Tcond_ptr)];
data = num2cell([zeros(1,2);num]);
data(1,:)={sprintf('Sa(T*=%g)',Tcond),'lambda Sa'};
c = uicontextmenu;
uimenu(c,'Label','Copy data','Callback',{@data2clipboard_uimenu,data});
uimenu(c,'Label','Undock','Callback'   ,{@figure2clipboard_uimenu,handles.ax1});
set(handles.ax1,'uicontextmenu',c);

num = [T,uhs,gmpe,cms,cmss1,cmss2];
data = num2cell([zeros(1,6);num]);
data(1,:) = {'T(s)','UHS','GMPE','CMS','CMS+SIG','CMS-SIG'};
c = uicontextmenu;
uimenu(c,'Label','Copy data','Callback',{@data2clipboard_uimenu,data});
uimenu(c,'Label','Undock','Callback'   ,{@figure2clipboard_uimenu,handles.ax2});
set(handles.ax2,'uicontextmenu',c);

c = uicontextmenu;
uimenu(c,'Label','Undock','Callback'   ,{@figure2clipboard_uimenu,handles.ax3});
set(handles.ax3,'uicontextmenu',c);

num = [T,rho];
data = num2cell([zeros(1,2);num]);
data(1,:)={'T(s)',sprintf('rho(T,T*)')};
c = uicontextmenu;
uimenu(c,'Label','Copy data','Callback',{@data2clipboard_uimenu,data});
uimenu(c,'Label','Undock','Callback'   ,{@figure2clipboard_uimenu,handles.ax4});
set(handles.ax4,'uicontextmenu',c);

format(cF);

end

function [handles,MRscen,rmin,rmax,Rcenter,Mcenter]=run_func(handles,sources,lambda2,deagg2)

deagg{1}    = vertcat(deagg2{1,1,1,:});
indsource   = zeros(0,2);
Nsources    = size(deagg2,4);
for i=1:Nsources
    dg         = deagg2{1,1,1,i};
    if ~isempty(dg)
        Nscen     = size(dg,1);
        NM        = length(sources(i).mscl(:,1));
        NR        = size(sources(i).aream,1);
        ptr2      = sort(repmat((1:NR)',NM,1));
        indsource = [indsource;[i*ones(Nscen,1),ptr2]];  %#ok<*AGROW>
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

Mbar = sum(deagg{1}(:,1).*deagg{1}(:,3))/sum(deagg{1}(:,3));
Rbar = sum(deagg{1}(:,2).*deagg{1}(:,3))/sum(deagg{1}(:,3));
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