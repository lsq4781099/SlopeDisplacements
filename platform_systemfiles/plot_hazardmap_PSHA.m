function plot_hazardmap_PSHA(fig,v,MapOptions,h,idx)

ax1=findall(fig,'Tag','ax1');
delete(findall(fig,'Tag','Colorbar'));
pall   = {'parula','autumn','bone','colorcube','cool','copper','flag','gray','hot','hsv','jet','lines','pink','prism','spring','summer','white','winter'};
fig.Colormap=feval(pall{MapOptions.map(2)});
ch=findall(ax1,'tag','gmap'); uistack(ch,'bottom');

nt = size(h.t,1);
for i=1:nt
    ch = findall(ax1,'tag',num2str(i));
    delete(ch);
end
ch=findall(ax1,'Tag','satext');delete(ch);

umax = -inf;
umin = inf;
label = h.id;
if ~isempty(h.shape)
    active = vertcat(h.shape.active);
    active = find(active);
end
ii     = 1;
text(NaN,NaN,'','parent',ax1,'Tag','satext');

% This piece of code prevents the patches to be generated for paths elements
t      = h.t;
ispath = regexp(t(:,1),'path');
for i=1:length(ispath)
    if isempty(ispath{i}),ispath{i}=0;end
end
ispath=cell2mat(ispath);
t(ispath==1,:)=[];
TriScatterd=cell(size(t,1),1);
for i=1:size(t,1)
    vptr=regexp(label,t{i});
    for j=1:length(vptr)
        if isempty(vptr{j})
            vptr{j}=0;
        end
    end
    vptr = find(cell2mat(vptr));
    x    = h.p(vptr,1);
    y    = h.p(vptr,2);
    
    if isempty(idx)
        u    = v(vptr);
    else
        u    = v(idx(vptr));
    end
    
    TriScatterd{i}=scatteredInterpolant(x,y,u,'linear','none');
    if regexp(t{i,1},'grid')
        conn = t{i,2};
        gps  = [x,y];
        in   = 1:size(gps,1);
    end
    
    if regexp(t{i,1},'shape')
        ind = active(ii);
        xl=h.shape(ind).Lon';
        yl=h.shape(ind).Lat';
        faces = h.shape(ind).faces;
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
        u = TriScatterd{i}(gps(:,1),gps(:,2));
        in = inpolygon(gps(:,1),gps(:,2),xl,yl);
        ii=ii+1;
    end
    shading = patch(...
        'parent',ax1,...
        'vertices',gps,...
        'faces',conn,...
        'facevertexcdata',u,...
        'facecol','interp',...
        'edgecol','none',...
        'linewidth',0.5,...
        'facealpha',0.7,...
        'Tag',num2str(i),...
        'ButtonDownFcn',{@site_click_PSHA,fig,h.p,MapOptions,TriScatterd,2},...
        'visible','on');
       
    uistack(shading,'bottom') % move map to bottom (so it doesn't hide previously drawn annotations)
    
    uin = u;
    uin = uin(in);
    umin = min(min(uin),umin);
    umax = max(max(uin),umax);
    
end

if umin==umax
    umin=umin-1e-3;
    umax=umax+1e-3;
end

if ~isempty(t)
    caxis([umin umax])
end

if ~isempty(t)
    ch=findall(fig,'tag','po_contours');ch.Enable='on';ch.Value=1;
    cb=colorbar('peer',ax1,'location','eastoutside','position',[0.94 0.16 0.02 0.65],'ylim',[umin,umax]);
    
    ch=findall(fig,'tag','IM_select');
    set(get(cb,'Title'),'String',ch.String{ch.Value})
    ax1.ButtonDownFcn={@clear_satxt;fig};
end

h = findall(ax1,'tag','gmap');
if ~isempty(h)
    h.ButtonDownFcn={@clear_satxt;fig};
    uistack(h,'bottom')
end

function clear_satxt(hObject,eventdata,fig) %#ok<*INUSL,*INUSD>
ch=findall(fig,'Tag','satext');delete(ch);

