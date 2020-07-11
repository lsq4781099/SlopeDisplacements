function[]=plotDEM(vertices,faces,value)

fig=figure;
fig.MenuBar='none';
ax = gca;
set(ax,'box','on','nextplot','add','fontsize',9)
rotate3d(ax);

patch(...
    'parent',ax,...
    'vertices',vertices,...
    'faces',faces,...
    'FaceVertexCData',value,...
    'FaceColor','interp',...
    'edgecolor','none');
xlabel(ax,'Lon','fontsize',9)
ylabel(ax,'Lat','fontsize',9)
H=colorbar;
set(get(H,'title'),'string','Z (m)')

