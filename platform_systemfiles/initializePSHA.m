
%% platform data structure
% general structure
handles.platformMode   = 1;
handles.pdffile        = [];

% site definitions
handles.h              = createObj('site');
handles.hc             = createObj('site');
handles.idx            = [];
handles.site_selection = [];
handles.site_colors    = zeros(0,3);

% analysis options
handles.sys            = [];
handles.opt            = createObj('opt');
handles.HazOptions     = createObj('hazoptions');

% psha parameters
handles.MRE            = [];
handles.MREPCE         = cell(0,0);

% deterministic scenarios
handles.source         = [];
handles.scenarios      = zeros(0,9);
handles.mulogIM        = [];
handles.L              = [];
handles.krate          = []; %k-mean cluster rate
handles.kY             = []; %k-mean cluster scenarios
handles.Y              = [];

%---initialization---------------------------------------------------
delete(findall(handles.fig,'type','legend'));
delete(findall(handles.fig,'tag','Colorbar'));
delete(findall(handles.fig,'type','line'));
delete(findall(handles.ax1,'type','patch'));
delete(findall(handles.ax1,'tag','gmap'));
delete(findall(handles.ax1,'type','text'));
handles.switchmode.CData=handles.form1;
xlabel(handles.ax1,'Lon','fontsize',8,'fontname','arial','tag','xlabel')
ylabel(handles.ax1,'Lat','fontsize',8,'fontname','arial','tag','ylabel')
set(handles.ax1,'nextplot','add','box','on','dataaspectratio',[1 1 1],'XGrid','off','YGrid','off');
handles.fig.Name='SeismicHazard - BETA';
handles.Boundary_check.Value=0;
handles.Layers_check.Value=0;
handles.site_menu.Value=1;
handles.site_menu.String={''};
handles.site_menu_psda.Value=1;
handles.site_menu_psda.String={''};
handles.IM_select.Value=1;
handles.IM_select.String={''};
handles.po_region.String='';
handles.po_region.Value=1;
handles.po_grid.Value=0;
handles.po_sources.Value=1;
handles.SourceLabels.Value=0;
handles.po_sourcemesh.Value=0;
handles.po_googleearth.Value=1;
handles.po_clusters.Value=0;
handles.po_region.Enable='off';
handles.po_region.Visible='off';
handles.PSHApannel.Visible='off';
handles.DSHApannel.Visible='off';
handles.engine.Enable='on';
handles.po_clusters.Enable='off';
handles.OpenRef.Visible='off';
handles.SeismicHazardTools.Enable='off';
handles.Edit.Enable='off';
handles.PSHAMenu.Enable='off';
handles.DSHAMenu.Enable='off';
handles.po_grid.Enable='on';
handles.po_sources.Enable='off';
handles.SourceLabels.Enable='off';
handles.po_sourcemesh.Enable='off';
handles.Boundary_check.Enable='off';
handles.Layers_check.Enable='off';
DATA=load('pshatoolbox_emptyGEmap.mat');
handles.ax1.XLim=DATA.XLIM;
handles.ax1.YLim=DATA.YLIM;
handles.ax1DefaultLimits=[DATA.XLIM;DATA.YLIM];
gmap=image(DATA.xx,DATA.yy,DATA.cc,'Parent',handles.ax1,'Tag','gmap','Visible','on');
handles.ax1.YDir='normal';
uistack(gmap,'bottom');
handles.ax1.Layer='top';
set(handles.ax2,'box','on','xscale','log','yscale','log')
xlabel(handles.ax2,'IM','fontsize',9)
ylabel(handles.ax2,'MeanRateofExceedance','fontsize',9)
plot_sites_PSHA;
setTickLabel(handles.ax1)
