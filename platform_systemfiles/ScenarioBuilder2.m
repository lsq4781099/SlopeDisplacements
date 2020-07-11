function varargout = ScenarioBuilder2(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @ScenarioBuilder2_OpeningFcn, ...
    'gui_OutputFcn',  @ScenarioBuilder2_OutputFcn, ...
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

function ScenarioBuilder2_OpeningFcn(hObject, ~, handles, varargin)
handles.output = [];
load All_Scenario_Buttons %#ok<LOAD>
handles.t.Data                  = cell(0,9);
handles.Undock.CData            = Undockbutton;
handles.ExitButton.CData        = ExitButtonbutton;
handles.plotTessel.CData        = plotTesselbutton;
handles.GridManager.CData       = GridManagerbutton;
handles.ShowNodes.CData         = ShowNodesbutton;
handles.Ruler.CData             = Rulerbutton;
handles.wiz1.CData              = invokeWizbutton;
handles.go_1.CData              = go_1;
handles.go_2.CData              = go_2;
handles.go_3.CData              = go_3;
handles.go_4.CData              = go_4;

handles.area       = patch('parent',handles.ax1,'faces',[],'vertices',[],'facecolor',[0.8 0.6 0.6],'facealpha',0.5,'tag','area','linewidth',1);
handles.edge6      = plot(handles.ax1,nan,nan,'-','color',[0 0 0],'tag','edge6','visible','on','linewidth',1);
%handles.RA         = patch('parent',handles.ax1,'faces',[],'vertices',[],'facecolor',[0.7 0.7 0.7],'facealpha',0.5,'tag','RA');
handles.meshR      = plot(handles.ax1,nan,nan,'+','markeredgecolor',[0.85 0.325 0.098],'tag','meshR','visible','on','markersize',6);
handles.hyp        = plot(handles.ax1,nan,nan,'ks','tag','area','visible','on','markersize',6,'markerfacecolor','r');
xlabel(handles.ax1,'X')
ylabel(handles.ax1,'Y')

% pre-proces input data
handles.sys        = varargin{1}; handles.branchmenu.String=compose('Branch %g',1:size(handles.sys.branch,1));
handles.opt        = varargin{2};
handles.h          = varargin{3};
handles.scenarios  = varargin{4};
handles.source     = varargin{5};
handles.pagenumber = 1;
handles.pages      = zeros(0,2);

% POPULATE TABLE WITH DATA
if ~isempty(handles.source)
    branch_ptr         = handles.source(1).branch;
    handles.branchmenu.Value =branch_ptr(1);
    
    Nscen              = size(handles.scenarios,1);
    handles.pages      = setDSHApages(Nscen);
    ind                = handles.pages(1,1):handles.pages(1,2);
    handles.t.Data     = num2cell(handles.scenarios(ind,1:9));
    DataR              = cell2mat(handles.t.Data(1,:));
    handles.pagenumber = 1;
    sSCEN_plotSource2(handles,handles.source(DataR(1)),DataR,handles.opt.ellipsoid);
    
    handles.go_1.Enable='on';
    handles.go_2.Enable='on';
    handles.go_3.Enable='on';
    handles.go_4.Enable='on';
    handles.pagetext.String = sprintf('1 of %g',size(handles.pages,1));
end
guidata(hObject, handles);
uiwait(handles.figure1);

function varargout = ScenarioBuilder2_OutputFcn(~,~, handles)
    
varargout{1} = handles.scenarios;
varargout{2} = handles.source;

delete(handles.figure1)

function figure1_CloseRequestFcn(hObject, ~, ~)  %#ok<*DEFNU>
if isequal(get(hObject,'waitstatus'),'waiting')
    uiresume(hObject);
else
    delete(hObject);
end
% delete(hObject);

function branchmenu_Callback(hObject, ~, handles)
brptr = handles.branchmenu.Value;
wipeDSHAmodel
handles.branchmenu.Value = brptr;
guidata(hObject,handles);

function wiz1_Callback(hObject, ~, handles)
bptr  = handles.branchmenu.Value;                           % branch pointer
gptr  = handles.sys.branch(bptr,1);                         % geometry pointer
slib  = handles.sys.labelG{gptr};
[sptr,tf] = listdlg('ListString',slib); % source pointer
if tf==0,return,end

% initializes Data
handles.source    = buildDSHAmodel2(handles.sys,handles.opt,handles.h,bptr,sptr);
handles.scenarios = sSCEN_list(handles.source);

% display magic
Nscen              = size(handles.scenarios,1);
handles.pages      = setDSHApages(Nscen);
handles.t.Data     = num2cell(handles.scenarios);
handles.pagenumber = 1;
DataR              = cell2mat(handles.t.Data(1,:));
sSCEN_plotSource2(handles,handles.source(DataR(1)),DataR,handles.opt.ellipsoid);

handles.go_1.Enable='on';
handles.go_2.Enable='on';
handles.go_3.Enable='on';
handles.go_4.Enable='on';
handles.pagetext.String = sprintf('1 of %g',size(handles.pages,1));

guidata(hObject,handles)

% ---------------- FILE MENU ---------------------------------------------
function File_Callback(~, ~, ~)

function Reset_Callback(hObject, ~, handles)
wipeDSHAmodel
guidata(hObject,handles);

function Exit_Callback(~,~, ~)
close(gcf)

function ExitButton_Callback(~,~, ~)
close(gcf)

% ------------ MANAGE SCENARIOS  -----------------------------------------
function ShowNodes_Callback(~, ~, handles)

switch handles.meshR.Visible
    case 'on'  , handles.meshR.Visible='off'; handles.voroni.Visible='off'; handles.voroni2.Visible='off';
    case 'off' , handles.meshR.Visible='on';  handles.voroni.Visible='on';  handles.voroni2.Visible='on';
end

function plotTessel_Callback(~, ~, handles)

% switch handles.voroni.Visible
%     case 'on'  , handles.voroni.Visible='off';handles.voroni2.Visible='off';
%     case 'off' , handles.voroni.Visible='on'; handles.voroni2.Visible='on';
% end

function GridManager_Callback(~, ~, handles)
xg = handles.ax1.XGrid;
switch xg
    case 'on' , handles.ax1.XGrid='off';handles.ax1.YGrid='off';
    case 'off', handles.ax1.XGrid='on';handles.ax1.YGrid='on';
end

% ------------- EDIT MENU ----------------------------------------------
function Edit_Callback(~, ~, ~)

% ------------- TABLE EDIT ------------------------------------------------
function Ruler_Callback(hObject, ~, handles)

ch1=findall(handles.ax1,'tag','patchselect');
ch2=findall(handles.ax1,'tag','patchtxt');
if isempty(ch1) && isempty(ch2)
    ch1=findall(handles.figure1,'Style','pushbutton','Enable','on'); set(ch1,'Enable','inactive');
    ch2=findall(handles.figure1,'type','uimenu','Enable','on'); set(ch2,'Enable','off');
    XYLIM1 = get(handles.ax1,{'xlim','ylim'});
    show_distanceECEF(handles.ax1,'line');
    XYLIM2 = get(handles.ax1,{'xlim','ylim'});
    set(handles.ax1,{'xlim','ylim'},XYLIM1);
    set(handles.ax1,{'xlim','ylim'},XYLIM2);
    set(ch1,'Enable','on')
    set(ch2,'Enable','on')
else
    delete(ch1);
    delete(ch2);
end
guidata(hObject, handles);

function Undock_Callback(hObject, eventdata, handles)
figure2clipboard_uimenu(hObject, eventdata,handles.ax1)

function t_CellSelectionCallback(hObject, eventdata, handles)
if isempty(eventdata.Indices)
    return
end
row      = eventdata.Indices(1);
DataR    = cell2mat(handles.t.Data(row,:));
sSCEN_plotSource2(handles,handles.source(DataR(1)),DataR,handles.opt.ellipsoid);
guidata(hObject,handles)

function go_1_Callback(hObject, ~, handles)
handles.pagenumber = 1;
ind                = handles.pages(1,1):handles.pages(1,2);
handles.t.Data     = num2cell(handles.scenarios(ind,1:8));
DataR              = cell2mat(handles.t.Data(1,:));
sSCEN_plotSource2(handles,handles.source(DataR(1)),DataR,handles.opt.ellipsoid);
handles.pagetext.String = sprintf('%g of %g',handles.pagenumber,size(handles.pages,1));
guidata(hObject,handles)

function go_2_Callback(hObject, ~, handles)
handles.pagenumber = max(handles.pagenumber-1,1);
ind                = handles.pages(handles.pagenumber,1):handles.pages(handles.pagenumber,2);
handles.t.Data     = num2cell(handles.scenarios(ind,1:8));
DataR              = cell2mat(handles.t.Data(1,:));
sSCEN_plotSource2(handles,handles.source(DataR(1)),DataR,handles.opt.ellipsoid);
handles.pagetext.String = sprintf('%g of %g',handles.pagenumber,size(handles.pages,1));
guidata(hObject,handles)

function go_3_Callback(hObject, ~, handles)

NumPages           = size(handles.pages,1);
handles.pagenumber = min(handles.pagenumber+1,NumPages);
ind                = handles.pages(handles.pagenumber,1):handles.pages(handles.pagenumber,2);
handles.t.Data     = num2cell(handles.scenarios(ind,1:8));
DataR              = cell2mat(handles.t.Data(1,:));
sSCEN_plotSource2(handles,handles.source(DataR(1)),DataR,handles.opt.ellipsoid);
handles.pagetext.String = sprintf('%g of %g',handles.pagenumber,size(handles.pages,1));
guidata(hObject,handles)

function go_4_Callback(hObject, ~, handles)

handles.pagenumber = size(handles.pages,1);
ind                = handles.pages(end,1):handles.pages(end,2);
handles.t.Data     = num2cell(handles.scenarios(ind,1:8));
DataR              = cell2mat(handles.t.Data(1,:));
sSCEN_plotSource2(handles,handles.source(DataR(1)),DataR,handles.opt.ellipsoid);
handles.pagetext.String = sprintf('%g of %g',handles.pagenumber,size(handles.pages,1));
guidata(hObject,handles)
