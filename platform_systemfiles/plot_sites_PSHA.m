delete(findall(handles.ax1,'tag','siteplot'));
delete(findall(handles.ax1,'tag','cluster'));
if isempty(handles.h.p)
    handles.po_sites.Enable='off';
    handles.po_clusters.Enable='off';
    handles.runMRE.Enable='off';
    return
else
    handles.po_sites.Enable='on';
    if strcmp(handles.opt.Clusters{1},'on')
        handles.po_clusters.Enable='on';
    else
        handles.po_clusters.Enable='off';
    end
end

c = repmat([0 0.4470 0.7410],size(handles.h.p,1),1);
switch handles.po_sites.Value
    case 0, scatter(handles.ax1,handles.h.p(:,1),handles.h.p(:,2),16,c,'filled','markeredgecolor','none','tag','siteplot','visible','off');
    case 1, scatter(handles.ax1,handles.h.p(:,1),handles.h.p(:,2),16,c,'filled','markeredgecolor','none','tag','siteplot','visible','on');
end

c = repmat([1 0 0],size(handles.hc.p,1),1);
switch and(strcmp(handles.opt.Clusters{1},'on'),handles.po_clusters.Value==1)
    case 0,scatter(handles.ax1,handles.hc.p(:,1),handles.hc.p(:,2),16,c,'filled','markeredgecolor','k','tag','cluster','visible','off');
    case 1,scatter(handles.ax1,handles.hc.p(:,1),handles.hc.p(:,2),16,c,'filled','markeredgecolor','k','tag','cluster','visible','on');
end

if isfield(handles,'sys')
    switch handles.sys.filename
        case 'USGS_NHSM_2008'
            handles.PSHAMenu.Enable='on';
            handles.DSHAMenu.Enable='off';
            handles.SeismicHazardTools.Enable='off';
            handles.runMRE.Enable='on';
        case 'USGS_NHSM_2014'
            handles.PSHAMenu.Enable='on';
            handles.DSHAMenu.Enable='off';
            handles.SeismicHazardTools.Enable='off';
            handles.runMRE.Enable='on';
        otherwise
            handles.PSHAMenu.Enable='on';
            handles.DSHAMenu.Enable='on';
            handles.SeismicHazardTools.Enable='on';
            handles.runMRE.Enable='on';
    end
else
    handles.Run.Enable='off';
end

if ~isempty(handles.h.id)
    handles.site_menu.String = handles.h.id;
    handles.site_menu_psda.String = handles.h.id;
    if isempty(handles.site_selection)
        handles.site_menu.Value  = 1;
        handles.site_menu_psda.Value  = 1;
    else
        handles.site_menu.Value       = handles.site_selection(1);
        handles.site_menu_psda.Value  = handles.site_selection(1);
    end
    ch=findall(handles.ax1,'tag','siteplot');
    uistack(ch, 'top');
end

