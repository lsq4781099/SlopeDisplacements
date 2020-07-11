function plotsiteselect(handles)

site  = cell2mat(handles.p.Data(:,2:end));
delete(findall(handles.ax1,'tag','siteplot'))
delete(findall(handles.figure2,'Tag','Colorbar'));
if isempty(site)
    return
end

Lon   = site(:,1)';
Lat   = site(:,2)';

val   = handles.layerpop.Value;
fld   = handles.layerpop.String{val};
CData = site(:,3+val);

scatter(handles.ax1,Lon,Lat,20,CData(:),'filled','markeredgecolor','k','tag','siteplot');
caxis('auto')
colormap(parula);
C=colorbar('peer',handles.ax1,'location','eastoutside','position',[0.94 0.16 0.02 0.65]);
C.Title.String=fld;
C.LimitsMode='auto';
