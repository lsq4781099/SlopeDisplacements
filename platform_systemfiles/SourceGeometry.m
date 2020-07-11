function varargout = SourceGeometry(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @SourceGeometry_OpeningFcn, ...
    'gui_OutputFcn',  @SourceGeometry_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end
if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end

function SourceGeometry_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>
handles.output    = hObject;
ihandles          = varargin{1};

icon = double(imread('exit.jpg'))/255; set(handles.Exit_button,   'CData',icon);
icon = double(imresize(imread('Ruler.jpg'),[16 16]))/255;     set(handles.setXYLimits,'CData',icon);

handles.opt       = ihandles.opt;
handles.sys       = ihandles.sys;
copyobj(findall(ihandles.ax1,'tag','gmap'), handles.ax1);
copyobj(findall(ihandles.ax1,'tag','shape1'), handles.ax1);
copyobj(findall(ihandles.ax1,'tag','shape2'), handles.ax1);

ch=findall(handles.ax1,'tag','gmap');   if ~isempty(ch),ch.Visible='off';end
ch=findall(handles.ax1,'tag','shape1'); if ~isempty(ch),ch.Visible='on';end
ch=findall(handles.ax1,'tag','shape2'); if ~isempty(ch),ch.Visible='off';end

handles.ax1.DataAspectRatio=ihandles.ax1.DataAspectRatio;
handles.ax1.XLim=ihandles.ax1.XLim;
handles.ax1.YLim=ihandles.ax1.YLim;

if isempty(handles.opt.ellipsoid.Code)
    xlabel(handles.ax1,'X (km)','fontsize',8,'fontname','arial');
    ylabel(handles.ax1,'Y (km)','fontsize',8,'fontname','arial')
    handles.tabla.ColumnName(1:2)={'Y(km)';'X(km)'};
else
    xlabel(handles.ax1,'Lon','fontsize',8,'fontname','arial');
    ylabel(handles.ax1,'Lat','fontsize',8,'fontname','arial')
    handles.tabla.ColumnName(1:3)={'Lat';'Lon';'Depth'};
end

Ngeom = size(handles.sys.Nsrc,2);
set(handles.po_region,'string',compose('Geometry %g',1:Ngeom)','value',1,'enable','on')
set(handles.source_list,'string',handles.sys.labelG{1},'value',1,'enable','on');
plot_geometry_PSHA(handles.ax1,handles.sys,handles.opt);
ngeom = length(handles.sys.src1);
for i=1:ngeom
    if i==handles.po_region.Value
        tag = sprintf('edge%g',i); ch  = findall(handles.ax1,'tag',tag); set(ch,'visible','on');
        tag = sprintf('mesh%g',i); ch  = findall(handles.ax1,'tag',tag); set(ch,'visible','on');
    else
        tag = sprintf('edge%g',i); ch  = findall(handles.ax1,'tag',tag); set(ch,'visible','off');
        tag = sprintf('mesh%g',i); ch  = findall(handles.ax1,'tag',tag); set(ch,'visible','off');
    end
end

plot_geometry_PSHA_single(handles);
rotate3d(handles.ax2);
xlabel(handles.ax2,'X(km)','fontsize',8)
ylabel(handles.ax2,'Y(km)','fontsize',8)
zlabel(handles.ax2,'Z(km)','fontsize',8)
handles.ax2.Color='none';
grid(handles.ax2,'off');
handles.ax1.Layer='top';
guidata(hObject, handles);
% uiwait(handles.figure1);

function varargout = SourceGeometry_OutputFcn(hObject, eventdata, handles)

varargout{1}=[];
% delete(handles.figure1)

function po_region_Callback(hObject, eventdata, handles) %#ok<*DEFNU>

ngeom = length(handles.sys.src1);
for i=1:ngeom
    if i==handles.po_region.Value
        tag = sprintf('edge%g',i); ch  = findall(handles.ax1,'tag',tag); set(ch,'visible','on');
        tag = sprintf('mesh%g',i); ch  = findall(handles.ax1,'tag',tag); set(ch,'visible','on');
    else
        tag = sprintf('edge%g',i); ch  = findall(handles.ax1,'tag',tag); set(ch,'visible','off');
        tag = sprintf('mesh%g',i); ch  = findall(handles.ax1,'tag',tag); set(ch,'visible','off');
    end
end

% update source list
handles.source_list.Value=1;
handles.source_list.String=handles.str{handles.po_region.Value};
plot_geometry_PSHA_single(handles);
guidata(hObject, handles);

function po_region_CreateFcn(hObject, eventdata, handles) %#ok<*INUSD>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function source_list_Callback(hObject, eventdata, handles)
plot_geometry_PSHA_single(handles);

function source_list_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Exit_button_Callback(hObject, eventdata, handles)
close(handles.figure1)

function Exit_Callback(hObject, eventdata, handles)
close(handles.figure1);

function tabla_CellSelectionCallback(hObject, eventdata, handles)
delete(findall(handles.ax1,'Tag','sourcepoint'));
delete(findall(handles.ax2,'Tag','sourcepoint'));
if isempty(eventdata.Indices)
    return
end
vertices = cell2mat(hObject.Data);
ind = eventdata.Indices(1);
plot(handles.ax1,vertices(ind,1),vertices(ind,2),'ko','markersize',5,'Tag','sourcepoint','markerfacecolor',[0.85 0.325 0.098]);

vertices = gps2xyz(vertices,handles.opt.ellipsoid);
plot3(handles.ax2,vertices(ind,1),vertices(ind,2),vertices(ind,3),'ko','markersize',5,'Tag','sourcepoint','markerfacecolor',[0.85 0.325 0.098]);
guidata(hObject, handles);

function tabla_CellEditCallback(hObject, eventdata, handles)

% store changes
pop   = get(handles.po_region,'value'); val   = get(handles.source_list,'value');
handles.model(pop).source(val).vertices=get(handles.tabla,'data');

ch=findall(handles.ax1,'Tag','sourcepoint');
vertices = get(hObject,'data');
ind = eventdata.Indices;
ch.XData=vertices(ind(1),2);
ch.YData=vertices(ind(1),1);

ch=findall(handles.ax2,'Tag','sourcepoint');
ch.XData=vertices(ind(1),2);
ch.YData=vertices(ind(1),1);
ch.ZData=vertices(ind(1),3);

guidata(hObject,handles)

function figure1_CloseRequestFcn(hObject, eventdata, handles)

% if isequal(get(hObject,'waitstatus'),'waiting')
%     uiresume(hObject);
% else
%     delete(hObject);
% end
delete(hObject);

function setXYLimits_Callback(hObject, eventdata, handles)

function File_Callback(hObject, eventdata, handles)

function checkbox2_Callback(hObject, eventdata, handles)

function po_grid_Callback(hObject, eventdata, handles)
val = hObject.Value;
switch val
    case 1, grid(handles.ax1,'on')
    case 0, grid(handles.ax1,'off')
end

function checkbox4_Callback(hObject, eventdata, handles)

function po_terrain_Callback(hObject, eventdata, handles)
ch=findall(handles.ax1,'tag','gmap');
if isempty(ch)
    return
end
switch hObject.Value
    case 0, ch.Visible='off';
    case 1, ch.Visible='on';
end

function po_boundaries_Callback(hObject, eventdata, handles)
ch = findall(handles.ax1,'tag','shape1');
val = hObject.Value;
switch val
    case 1, set(ch,'visible','on');
    case 0, set(ch,'visible','off')
end
guidata(hObject,handles);

function po_sourcemesh_Callback(hObject, eventdata, handles)

i   = handles.po_region.Value;
tag_mesh = sprintf('mesh%g',handles.po_region.Value);
switch hObject.Value
    case 1
        ch=findall(handles.ax1,'tag',tag_mesh); set(ch,'visible','on');
    case 0
        ch=findall(handles.ax1,'tag',tag_mesh); set(ch,'visible','off');
end
guidata(hObject,handles);

function po_layers_Callback(hObject, eventdata, handles)

function po_sourcelabels_Callback(hObject, eventdata, handles)

val = handles.po_region.Value;
tag = sprintf('sourcetag%g',val);
ch  = findall(handles.ax1,'tag',tag);
switch hObject.Value
    case 0, set(ch,'visible','off')
    case 1, set(ch,'visible','on')
end

function Undock_Callback(hObject, eventdata, handles)

function SincleSourceUndock_Callback(hObject, eventdata, handles)
figure2clipboard_uimenu(hObject, eventdata,handles.ax2)

function Allsources_Callback(hObject, eventdata, handles)
figure2clipboard_uimenu(hObject, eventdata,handles.ax1)
