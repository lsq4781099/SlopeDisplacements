
handles.site_selection = 1:size(handles.h.p,1);
handles.fig.Name=[handles.sys.filename,' - Seismic Hazard'];
delete(findall(handles.ax2,'type','line'));

% ---  Updates Period Selection ---
handles.IM_select.Value  = 1;
handles.IM_select.String = IM2str(handles.opt.IM);

% Updates Po-Region
ngeom =size(handles.sys.Nsrc,2);
if ngeom>1
    handles.po_region.Enable='on';
    handles.po_region.Visible='on';
    handles.po_region.String=compose('Geometry %i',1:ngeom);
end

plot_geometry_PSHA(handles.ax1,handles.sys,handles.opt);
if handles.po_sources.Value
    tag = sprintf('edge%g',handles.po_region.Value); ch=findall(handles.ax1,'tag',tag); set(ch,'visible','on');
end
if handles.po_sourcemesh.Value
    tag = sprintf('mesh%g',handles.po_region.Value); ch=findall(handles.ax1,'tag',tag); set(ch,'visible','on');
end

xlabel(handles.ax2,addIMunits(handles.IM_select.String{1}),'fontsize',8);
handles.po_sources.Enable        = 'on';
handles.po_sourcemesh.Enable     = 'on';
handles.SourceLabels.Enable      = 'on';
handles.InspectInputFile.Enable  = 'on';
handles.Edit.Enable              = 'on';
handles.Boundary_check.Enable    = 'on'; handles.Boundary_check.Value = 1;
handles.Layers_check.Enable      = 'on'; handles.Layers_check.Value   = 0;
default_maps(handles.fig,handles.opt);

if isempty(handles.opt.Image)
    rat = diff(handles.ax1.YLim)/diff(handles.ax1.XLim);
    axis(handles.ax1,'equal');
    axis(handles.ax1,'auto');
    XL = handles.ax1.XLim;
    YC = mean(handles.ax1.YLim);
    handles.ax1.YLim=YC+rat*diff(XL)*[-1 1];
end

if isfield(handles,'h') && ~isempty(handles.h.id)
    handles.PSHAMenu.Enable='on';
    handles.DSHAMenu.Enable='on';
    handles.SeismicHazardTools.Enable='on';
    handles.runMRE.Enable='on';
end

if strcmp(handles.opt.Clusters{1},'on')
    handles.po_clusters.Enable='on';
end
[handles.idx,handles.hc] = compute_clusters(handles.opt,handles.h);
Nsites = size(handles.h.p,1);
handles.site_colors = repmat([0 0 1],Nsites,1);


handles.ax1DefaultLimits=[handles.ax1.XLim;handles.ax1.YLim];