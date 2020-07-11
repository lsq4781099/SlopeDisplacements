function varargout = LIBS(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @LIBS_OpeningFcn, ...
    'gui_OutputFcn',  @LIBS_OutputFcn, ...
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

function LIBS_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>

load libsbuttons c1 c2 c3 c4 c5 c6 c7 c8 c9
handles.CDataOpen                 = c1;
handles.CDataClosed               = c2;
handles.toggle1.CData             = c3;
handles.runREG.CData              = c4;
handles.ax2Limits.CData           = c5;
handles.REG_DisplayOptions.CData  = c6;
handles.treebutton.CData          = c7;
handles.deleteButton.CData        = c8;
handles.undock.CData              = c9;
xlabel(handles.ax1,'Settlement (mm)','fontsize',10)
ylabel(handles.ax1,'Mean Rate of Exceedance','fontsize',10)

%% Retrieve data from SeismicHazard
[handles.REG_Display] = defaultPSDA_plotoptions;
handles.haz   = [];
handles.haz2  = [];
handles       = wipeLIBS(handles);
guidata(hObject, handles);
% uiwait(handles.fig);

function runREG_Callback(hObject, eventdata, handles)
if ~isfield(handles,'sys'), return;end
handles.haz = haz_LIBS(handles);
handles     = runLIBS_regular(handles);
plot_LIBS_regular(handles);
handles.deleteButton.CData  = handles.CDataClosed;
guidata(hObject,handles)

function varargout = LIBS_OutputFcn(hObject, eventdata, handles)
varargout{1} = [];

% ----------------  FILE MENU ---------------------------------------------
function File_Callback(hObject, eventdata, handles) %#ok<*INUSD,*DEFNU>

function ImportSeismicHazard_Callback(hObject, eventdata, handles)

if isfield(handles,'defaultpath_others')
    [filename,pathname,FILTERINDEX]=uigetfile(fullfile(handles.defaultpath_others,'*.txt'),'select file');
else
    [filename,pathname,FILTERINDEX]=uigetfile(fullfile(pwd,'*.txt'),'select file');
end
if FILTERINDEX==0,return;end

handles.defaultpath_others=pathname;
handles = wipeLIBS(handles);
handles = importLIBS(handles,pathname,filename);

handles.fig.Name = sprintf('%s - Liquefaction Induced Building Settlement - LIBS',filename);
guidata(hObject,handles)

function Exit_Callback(hObject, eventdata, handles)
close(handles.fig)

function fig_CloseRequestFcn(hObject, eventdata, handles)
delete(hObject);

% ----------------  EDIT MENU ---------------------------------------------
function Edit_Callback(hObject, eventdata, handles)

function SModelExplorer_Callback(hObject, eventdata, handles)
SettleTiltExplorer

% ----------------  ANALYSIS PARAMETERS MENU ---------------------------------------------
function Analysis_Callback(hObject, eventdata, handles)

function deleteButton_Callback(hObject, eventdata, handles)
delete(findall(handles.ax1,'type','line'));
handles.deleteButton.CData  = handles.CDataOpen;
guidata(hObject,handles)

% ----------------  POP MENUS ---------------------------------------------

function pop_source_Callback(hObject, eventdata, handles)

function pop_source_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function pop_site_Callback(hObject, eventdata, handles)
plot_LIBS_regular(handles);
guidata(hObject,handles)

function pop_site_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% ----------------  PLOTTING  ---------------------------------------------

function mouseClick(hObject,eventdata,handles)

x        = getAbsCoords(handles.ax1);
d_ptr    = interp1(handles.d,(1:length(handles.d))',x,'nearest','extrap');
d        = handles.d(d_ptr);
ch       = findall(handles.ax1,'tag','line');
ch.XData = d*[1 1];
ch.YData = handles.ax1.YLim;

function [x, y] = getAbsCoords(h_ax)
crd = get(h_ax, 'CurrentPoint');
x = crd(2,1);
y = crd(2,2);

function tableREG_CellSelectionCallback(hObject, eventdata, handles)
guidata(hObject,handles)

function REG_DisplayOptions_Callback(hObject, eventdata, handles)
Nbranches = size(handles.tableREG.Data,1);
str       = compose('Branch %g',(1:Nbranches)');
handles.REG_Display = Display_Options_REG(str,handles.REG_Display);
plot_LIBS_regular(handles);
guidata(hObject,handles)

function treebutton_Callback(hObject, eventdata, handles)

if ~isfield(handles,'sys')
    return
end

[handles.T1,handles.T2,handles.T3,handles.optlib]=LIBSLogicTree(...
    handles.setLIB,handles.T1,handles.T2,handles.T3,handles.optlib,handles.sys.weight(:,5));

[handles.tableREG.Data,handles.IJK]=main_settle(handles.T1,handles.T2,handles.T3);
guidata(hObject,handles)

function ax2Limits_Callback(hObject, eventdata, handles)
handles.ax1=ax2Control(handles.ax1);

function OptionsPSDA_Callback(hObject, eventdata, handles)
% if isfield(handles,'paramPSDA')
%     handles.paramPSDA=PSDA2Parameters(handles.paramPSDA);
%     guidata(hObject,handles)
% end

function ClearModel_Callback(hObject, eventdata, handles)

handles=wipeLIBS(handles);
guidata(hObject,handles)

function handles=wipeLIBS(handles)

if isfield(handles,'sys')
    handles=rmfield(handles,{'sys','opt','h'});
end
% [handles.REG_Display]   = defaultSETTLE_plotoptions;
handles.haz               = [];
handles.fig.Name          = 'Liquefaction Induced Building Settlement';
handles.tableREG.Data     = cell(0,4);
delete(findall(handles.ax1,'type','line'))

function Tools_Callback(hObject, eventdata, handles)

function toggle1_Callback(hObject, eventdata, handles)

ch = findall(handles.fig,'type','legend');
if ~isempty(ch)
    switch ch.Visible
        case 'on',  ch.Visible='off';
        case 'off', ch.Visible='on';
    end
end
guidata(hObject,handles);

function OpenDrivingFile_Callback(hObject, eventdata, handles)
if ispc
    winopen(handles.sys.filename)
end

function undock_Callback(hObject, eventdata, handles)
figure2clipboard_uimenu(hObject, eventdata,handles.ax1)
