function plot_geometry_PSHA_single(handles)

pop   = get(handles.po_region,'value');
val   = get(handles.source_list,'value');
numG  = handles.sys.numG{pop}(val,:);
sourcetype = numG(1);
localindex = numG(2);


delete(findall(handles.ax1,'tag','selectedobject'));
chd =findall(handles.ax2,'tag','selectedobject');
meshcolor  = [0.6 0.6 1];
ellip      = handles.opt.ellipsoid;

source     = handles.sys.(sprintf('src%g',sourcetype));
vptr       = source.vptr(localindex,1):source.vptr(localindex,2);

if sourcetype~=1
    mptr       = source.mptr(localindex,1):source.mptr(localindex,2);
    
end
vert       = source.vert(vptr,:);
viewpoint  = mean(vert,1);
viewpoint(3)=100;
viewpoint  = gps2xyz(viewpoint,ellip);
handles.tabla.Data=num2cell(vert);


switch sourcetype
    case 1
        vert2 = gps2xyz(vert,handles.opt.ellipsoid);
        plot3(handles.ax2,vert2(1) ,vert2(2) ,vert2(3),'ks','markerfacecolor','k','tag','selectedobject')
    case 2
        plot (handles.ax1,vert([1:end,1],1) ,vert([1:end,1],2),'k-','tag','selectedobject','linewidth',1)
        conn    = source.conn(mptr,:);
        patch(handles.ax2,'vertices',source.xyzm,'faces',conn,'marker','s','markerfacecolor','k','facecolor','none','edgecolor','k','tag','selectedobject');
    case {3,4,5,6}
        plot (handles.ax1,vert([1:end,1],1) ,vert([1:end,1],2),'k-','tag','selectedobject','linewidth',1)
        vert2   = gps2xyz(vert,handles.opt.ellipsoid);
        conn    = source.conn(mptr,:);
        patch(handles.ax2,'vertices',source.xyzm,'faces',conn,'facecolor','w','edgecolor',meshcolor,'tag','selectedobject');
        plot3(handles.ax2,vert2([1:end,1],1),vert2([1:end,1],2),vert2([1:end,1],3),'k-','tag','selectedobject','linewidth',1.5);
    case 7
        plot (handles.ax1,vert([1:end,1],1) ,vert([1:end,1],2),'k-','tag','selectedobject','linewidth',1)
        vert2   = gps2xyz(vert,handles.opt.ellipsoid);
        
        p = source.xyzm;
        t = source.conn(mptr,:);
        t2 = [t(:,[1 2 3]);t(:,[1 4 3]);t(:,[2 3 4]);t(:,[1 2 4])];
        t2 = unique(sort(t2,2),'rows');
        patch(handles.ax2,'vertices',p,'faces',t2,'facecolor','w','edgecolor',meshcolor,'tag','selectedobject');
        plot3(handles.ax2,vert2([1:end,1],1),vert2([1:end,1],2),vert2([1:end,1],3),'k-','tag','selectedobject','linewidth',1.5);
        
end

delete(chd)

if ~isempty(handles.opt.ellipsoid.Code)
    set(handles.ax2, 'CameraPosition', viewpoint);
end