function default_maps(fig,opt)
ax1            = findall(fig,'tag','ax1');
Boundary_check = findall(fig,'tag','Boundary_check');
Layers_check   = findall(fig,'tag','Layers_check');

delete(findall(ax1,'tag','gmap'));
delete(findall(ax1,'tag','shape1'));
delete(findall(ax1,'tag','shape2'));

if isempty(opt.Image) 
    opt.Image = 'pshatoolbox_emptyGEmap.mat';
end
DATA         = load(opt.Image);
ax1.XLim     = DATA.XLIM;
ax1.YLim     = DATA.YLIM;
gmap         = image(DATA.xx,DATA.yy,DATA.cc, 'Parent', ax1);
gmap.Tag     = 'gmap';
gmap.Visible = 'on';
ax1.YDir      = 'normal';
uistack(gmap,'bottom');

ax1.ColorOrderIndex=1;
% shape1
if ~isempty(opt.Boundary)
    data = shaperead(opt.Boundary, 'UseGeoCoords', true);
    Boundary_check.Enable='on';
    switch Boundary_check.Value
        case 0,vis = 'off';
        case 1,vis = 'on';
    end
    plot(ax1,horzcat(data.Lon),horzcat(data.Lat),'tag','shape1','visible',vis);
end

ax1.Layer='top';
setTickLabel(ax1);

