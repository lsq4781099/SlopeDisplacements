
pd = makedist('normal');
sigma = handles.opt.Sigma;
if ~isempty(sigma) && strcmpi(sigma{1},'truncate')
    pd = truncate(pd,-inf,sigma{2});
elseif ~isempty(sigma) && strcmpi(sigma{1},'overwrite')
    pd.sigma=sigma{2};
end

uscen    = unique(handles.scenarios(:,1:2),'rows','stable');
scen     = handles.scenarios(handles.pop_scenario.Value,:);

if size(handles.L,3)==1
    ind2=1;
else
    [~,ind2] = intersect(uscen,scen(1:2),'rows');
end
source   = handles.source(scen(1));
mu       = permute(handles.mulogIM(handles.pop_scenario.Value,:,:),[3 2 1]);
Ls       = handles.L(:,:,ind2);

ax1 = findall(handles.fig,'tag','ax1');
ax2 = findall(handles.fig,'tag','ax2');

delete(findall(ax1,'tag','scenario'));

ch2 = findall(handles.fig,'tag','po_sites');
switch ch2.Value
    case 0, VIS1='off';
    case 1, VIS1='on';
end

ch3 = findall(handles.fig,'tag','po_contours');
switch ch3.Value
    case 0, VIS2='off';
    case 1, VIS2='on';
end

switch ch2.Value+ch3.Value
    case 0,    VIS3='off';
    otherwise, VIS3='on';
end

ch4 = findall(handles.fig,'tag','pop_field');
IM_ptr = ch4.Value;
if iscell(ch4.String)
    xlabel(ax2,ch4.String{IM_ptr},'fontsize',9)
else
    xlabel(ax2,ch4.String,'fontsize',9)
end
ylabel(ax2,'Mean Rate of Exceedance','fontsize',9)

% create scenario
Nsites  = length(handles.h.id);
II      = (1:Nsites)+Nsites*(IM_ptr-1);
mulogIM = mu(:,IM_ptr);
Nz      = size(Ls,1);

ch5 = findall(handles.fig,'tag','PSDA_display_mode');
switch ch5.Value
    case 1
        Y = exp(mulogIM);
    case 2
        Y = normrnd(0,1,Nsites,1);
        Z = random(pd,[Nsites,1]);
    case 3
        %Z = normrnd(0,1,Nz,1);
        Z = random(pd,[Nz,1]);
        Y = exp(mulogIM+Ls(II,:)*Z);
end
delete(findall(handles.fig,'Tag','Colorbar'));
pallet = 'parula';
set(handles.fig,'colormap',feval(pallet));

% This piece of code prevents the patches to be generated for paths elements
t      = handles.h.t;
ispath = regexp(t(:,1),'path');
for i=1:length(ispath)
    if isempty(ispath{i}),ispath{i}=0;end
end
ispath=cell2mat(ispath);
t(ispath==1,:)=[];

delete(findall(ax1,'tag','siteplot'));
delete(findall(ax1,'Tag','satext'));
XY      = handles.h.p(:,1:2);

if size(t,1)==0
    scatter(ax1,XY(:,1),XY(:,2),[],Y(:),'filled','markeredgecolor','k','tag','siteplot','visible',VIS1);
    cb=colorbar('peer',ax1,'location','eastoutside','position',[0.94 0.16 0.02 0.65]);
    CYLIM=sort([0.9*min(Y),1.1*max(Y)]);
    cb.YLim = CYLIM;
    caxis(ax1,CYLIM);
    cb.Visible=VIS3;
    return
else
    ch3.Enable='on';
end
% ------------------------------------------------------------------------

ch=findall(ax1,'tag','gmap'); uistack(ch,'bottom');
nt = size(handles.h.t,1);
for i=1:nt
    ch = findall(ax1,'tag',num2str(i));
    delete(ch);
end

ch=findall(ax1,'Tag','satext');delete(ch);
umax  = -inf;
umin  = inf;
label = handles.h.id;
p     = handles.h.p(:,[1,2]);

if ~isempty(handles.h.shape)
    active = vertcat(handles.h.shape.active);
    active = find(active);
end
ii     = 1;
text(NaN,NaN,'','parent',ax1,'Tag','satext');
for i=1:size(t,1)
    vptr=regexp(label,t{i});
    for j=1:length(vptr)
        if isempty(vptr{j})
            vptr{j}=0;
        end
    end
    vptr = find(cell2mat(vptr));
    x    = p(vptr,1);
    y    = p(vptr,2);
    u    = Y(vptr);
    F    = scatteredInterpolant(x,y,u,'linear','none');
    
    if regexp(t{i,1},'grid')
        conn = t{i,2};
        gps  = [x,y];
        in   = 1:size(gps,1);
    end
    
    if regexp(t{i,1},'shape')
        ind = active(ii);
        xl=handles.h.shape(ind).Lon';
        yl=handles.h.shape(ind).Lat';
        faces = handles.h.shape(ind).faces;
        pv =[xl,yl];
        gps  = zeros(0,2);
        conn = zeros(0,3);
        offset=0;
        for count=1:size(faces,1)
            indface = faces(count,:);
            indface(isnan(indface))=[];
            pvface = pv(indface,:);
            [gps_i,conn_i]=triangulate_maps(pvface,[x,y]);
            gps = [gps;gps_i]; %#ok<AGROW>
            conn = [conn;conn_i+offset]; %#ok<AGROW>
            offset=size(gps,1);
        end
        u = F(gps(:,1),gps(:,2));
        in = inpolygon(gps(:,1),gps(:,2),xl,yl);
        ii=ii+1;
    end
    
    
    
    shad=patch(...
        'parent',ax1,...
        'vertices',gps,...
        'faces',conn,...
        'facevertexcdata',u,...
        'facecol','interp',...
        'edgecol','none',...
        'linewidth',0.5,...
        'facealpha',0.7,...
        'Tag',num2str(i),...
        'visible',VIS2);
    uistack(shad,'bottom') % move map to bottom (so it doesn't hide previously drawn annotations)
    
    uin = u;
    uin = uin(in);
    umin = min(min(uin),umin);
    umax = max(max(uin),umax);
    
end

caxis([umin umax])

switch ch2.Value
    case 0, VIS='off';
    case 1, VIS='on';
end

%scatter(ax1,XY(:,1),XY(:,2),[],Y,'filled','markeredgecolor','k','tag','siteplot','visible',VIS1,'ButtonDownFcn',{@site_click_PSDA;handles;1},'visible',VIS);
scatter(ax1,XY(:,1),XY(:,2),[],Y(:),'filled','markeredgecolor','k','tag','siteplot','visible',VIS1,'visible',VIS);
cb=colorbar('peer',ax1,'location','eastoutside','position',[0.94 0.16 0.02 0.65],'YLim',[umin,umax]);
cb.Visible=VIS3;
%ax1.ButtonDownFcn={@clear_satxt,handles};

% restores map
gmap = findall(ax1,'tag','gmap');
if ~isempty(gmap)
    %handles.h.ButtonDownFcn={@clear_satxt;handles};
    uistack(gmap,'bottom') % move map to bottom (so it doesn't hide previously drawn annotations)
end

% function clear_satxt(hObject,eventdata,handles) %#ok<*INUSL,*INUSD>
% findall(ax1,'Tag','satext');

