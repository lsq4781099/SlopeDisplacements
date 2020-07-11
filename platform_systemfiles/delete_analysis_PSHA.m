handles.MRE       = [];
handles.MREPCE    = cell(0,0);
handles.Y         = [];
handles.source    = [];
handles.scenarios = zeros(0,9);
delete(findall(handles.fig,'tag','Colorbar'));
delete(findall(handles.ax2,'Type','line'));
delete(findall(handles.fig,'tag','Colorbar'));
delete(findall(handles.fig,'type','legend'));
delete(findall(handles.ax1,'Tag','scenario'));
delete(findall(handles.ax1,'Tag','satext'));
nt = size(handles.h.t,1);
for i=1:nt
    ch = findall(handles.ax1,'tag',num2str(i));
    delete(ch);
end

Nsites = size(handles.h.p,1);
if Nsites>0
    handles.site_colors = repmat([0 0 1],Nsites,1);
end