function sSCEN_plotSource2(handles,source,Data,ellipsoid) %#ok<*INUSD>

X       = Data(3);
Y       = Data(4);
Z       = Data(5);

ind     = sum(handles.scenarios(:,1:2)-Data(1,1:2),2)==0;
xx      = handles.scenarios(ind,3);
yy      = handles.scenarios(ind,4);
zz      = handles.scenarios(ind,5);


xlabel(handles.ax1,'X')
ylabel(handles.ax1,'Y')
vertices = source.vert;
Nnodes   = size(vertices,1);
gps      = xyz2gps([X,Y,Z],ellipsoid);
gpsm     = xyz2gps([xx,yy,zz],ellipsoid);


if source.obj==6
    gpsmA = xyz2gps(source.xyzm,ellipsoid);
    set(handles.area   ,'faces',source.conn(:,1:3),'vertices',gpsmA,'edgecolor','none');
    set(handles.edge6  ,'xdata',source.vert(:,1),'ydata',source.vert(:,2));
    set(handles.hyp    ,'XData',gps(1)   ,'YData',gps(2));
   set(handles.meshR  ,'XData',gpsm(:,1),'YData',gpsm(:,2));
else
    set(handles.area   ,'faces',1:Nnodes ,'vertices',vertices);
    set(handles.hyp    ,'XData',gps(1)   ,'YData',gps(2));
    set(handles.meshR  ,'XData',gpsm(:,1),'YData',gpsm(:,2));
    set(handles.edge6  ,'xdata',[],'ydata',[]);
end
