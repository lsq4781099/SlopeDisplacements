function varargout = SettleTiltExplorer(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @SettleTiltExplorer_OpeningFcn, ...
    'gui_OutputFcn',  @SettleTiltExplorer_OutputFcn, ...
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

function SettleTiltExplorer_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>
handles.output = hObject;
load icons_SettleTiltExplorer.mat c
handles.DiscretizeMButtonbutton.CData = c{1};
handles.Book.CData    = c{2};
handles.undock1.CData = c{3};
handles.undock2.CData = c{3};

methods = pshatoolbox_methods(6);
handles.func = {methods.func};
handles.STpop.String = {methods.label};
handles.STpop.Value  = 1;

% Write xlabel and ylabel
handles.plotmodeax2='ccdf';
handles.ax1XLABEL=xlabel(handles.ax1,'CAV (cm/s)');
handles.ax1YLABEL=ylabel(handles.ax1,'S (mm)');
handles.ax2XLABEL=xlabel(handles.ax2,'s (mm)');
handles.ax2YLABEL=ylabel(handles.ax2,'ccdf');

% fills Mass and contact pressure
handles=popSettleModel(handles,'Default_LIQ.txt');
plotSettle(handles)

guidata(hObject, handles);
% uiwait(handles.figure1);

function varargout = SettleTiltExplorer_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;

function handles=popSettleModel(handles,filename)

if nargin==2
    param=loadsiteLIBS(filename);
    handles.param = param;
    handles.B    .String = sprintf('%g',param.B   );
    handles.L    .String = sprintf('%g',param.L   );
    handles.Df   .String = sprintf('%g',param.Df  );
    handles.Q    .String = sprintf('%g',param.Q(1));
    [~,handles.bType.Value]= intersect(handles.bType.String,param.type);
    handles.wt   .String = sprintf('%g',param.wt  );
    
    
    [~,handles.LPC.Value]  = intersect(handles.LPC.String,param.LPC);
    handles.th1  .String = sprintf('%g',param.th1 );
    handles.th2  .String = sprintf('%g',param.th2 );
    handles.N1   .String = sprintf('%g',param.N1  );
    handles.N2   .String = sprintf('%g',param.N2  );
    [~,handles.meth.Value]  = intersect(handles.meth.String,param.meth);
    
    switch param.meth
        case 'SPT', handles.table1.Data = num2cell([param.N160(:),param.thick(:),param.d2mat(:)]);
        case 'CPT', handles.table1.Data = num2cell([param.qc1N(:),param.thick(:),param.d2mat(:)]);
    end
    
    handles.table2.Data = num2cell(param.CPT);
end

handles.text1.Visible='off'; handles.e1.Visible='off';
handles.text2.Visible='off'; handles.e2.Visible='off';
handles.text3.Visible='off'; handles.e3.Visible='off';
handles.text4.Visible='off'; handles.e4.Visible='off';
handles.pan1.Visible = 'off';
handles.pan2.Visible = 'off';

switch handles.STpop.Value
    case 1
        handles.pan1.Visible  = 'on';
        handles.text1.Visible = 'on';
        handles.e1.Visible    = 'on';
        
        handles.text1.String  = 'CAV (cm/s)';
        handles.e1.String     = 200;
    case 2
        handles.pan2.Visible  = 'on';
        handles.text1.Visible = 'on';
        handles.text2.Visible = 'on';
        handles.text3.Visible = 'on';
        handles.text4.Visible = 'on';
        handles.e1.Visible    = 'on';
        handles.e2.Visible    = 'on';
        handles.e3.Visible    = 'on';
        handles.e4.Visible    = 'on';
        
        handles.text1.String  = 'Mag';
        handles.text2.String  = 'PGA';
        handles.text3.String  = 'Sa1 (g)';
        handles.text4.String  = 'CAVdp (cm/s)';
        handles.e1.String     = 7.5;
        handles.e2.String     = 0.45;
        handles.e3.String     = 0.90;
        handles.e4.String     = 200;
    case 3
        handles.pan1.Visible  = 'on';
        handles.text1.Visible = 'on';
        handles.e1.Visible    = 'on';
        
        handles.text1.String  = 'CAV (cm/s)';
        handles.e1.String     = 200;
    case 4
        handles.pan1.Visible  = 'on';
        handles.text1.Visible = 'on';
        handles.text2.Visible = 'on';
        handles.e1.Visible    = 'on';
        handles.e2.Visible    = 'on';
        
        handles.text1.String  = 'CAV (cm/s)';
        handles.text2.String  = 'VGI (cm/s)';
        handles.e1.String     = 200;
        handles.e2.String     = 200;
end

function []=plotSettle(handles)

cav0  = logsp(1,2000,150); % CAV (cm/s)
param.B=str2double(handles.B.String);
param.L=str2double(handles.L.String);
param.Df=str2double(handles.Df.String);
param.Q=str2double(handles.Q.String);
param.type=handles.bType.String{handles.bType.Value};
param.wt=str2double(handles.wt.String);
param.LPC=handles.LPC.String{handles.LPC.Value};
param.th1=str2double(handles.th1.String);
param.th2=str2double(handles.th2.String);
param.N1=str2double(handles.N1.String);
param.N2=str2double(handles.N2.String);
param.meth='SPT';%handles.meth.String{handles.meth.Value};
param.N160=cell2mat(handles.table1.Data(:,1));
param.qc1N=cell2mat(handles.table1.Data(:,1));
param.thick=cell2mat(handles.table1.Data(:,2));
param.d2mat=cell2mat(handles.table1.Data(:,3));
param.CPT=handles.param.CPT;

switch handles.STpop.Value
    case 1
        s     = logsp(1,500,200);
        cav         = str2double(handles.e1.String);
        [lny0,sig0] = set_B18(param,cav0);
        [mu,sig]    = set_B18(param,cav);
        handles.ax1XLABEL.String='CAV (m/s)';
        handles.ax1YLABEL.String='S (mm)';
        handles.ax2XLABEL.String='s (mm)';
    case 2
        s     = logsp(1,500,200);
        M     = str2double(handles.e1.String);
        PGA   = str2double(handles.e2.String); 
        Sa1   = str2double(handles.e3.String);
        CAV   = str2double(handles.e4.String); 
        [LBS,HL] = calc_LBS_FS (param.CPT, param.wt, param.Df,M, PGA);
        handles.LBS.String=sprintf('%4.3g',LBS);
        handles.HL.String =sprintf('%3.2g',HL);
        
        param.LBS=LBS;
        param.HL =HL;
        [lny0,sig0] = set_I17(param,Sa1,cav0);
        [mu,sig]    = set_I17(param,Sa1,CAV);
        handles.ax1XLABEL.String='CAV (m/s)';
        handles.ax1YLABEL.String='S (mm)';
        handles.ax2XLABEL.String='s (mm)';
    case 3
        s           = logsp(0.01,10,200);
        cav         = str2double(handles.e1.String);
        [lny0,sig0] = til_B18e(param,cav0);
        [mu,sig]    = til_B18e(param,cav);
        handles.ax1XLABEL.String='CAV (m/s)';
        handles.ax1YLABEL.String='\theta (°)';
        handles.ax2XLABEL.String='\theta (°)';
    case 4
        s     = logsp(0.01,10,200);
        cav         = str2double(handles.e1.String);
        vgi0        = str2double(handles.e2.String);
        [lny0,sig0] = til_B18s(param,cav0,vgi0);
        [mu,sig]    = til_B18s(param,cav,vgi0);
        handles.ax1XLABEL.String='CAV (m/s)';
        handles.ax1YLABEL.String='\theta (°)';
        handles.ax2XLABEL.String='\theta (°)';
end

switch handles.plotmodeax2
    case 'pdf'  , Y   =   normpdf(log(s),mu,sig);  handles.ax2YLABEL.String='pdf';
    case 'cdf'  , Y   =   normcdf(log(s),mu,sig);  handles.ax2YLABEL.String='cdf';
    case 'ccdf' , Y   = 1-normcdf(log(s),mu,sig);  handles.ax2YLABEL.String='ccdf';
end

delete(findall(handles.ax1,'tag','curves'))
handles.ax1.ColorOrderIndex=1; plot(handles.ax1,cav0,exp(lny0)     ,'','tag','curves','DisplayName','\mu')
handles.ax1.ColorOrderIndex=1; plot(handles.ax1,cav0,exp(lny0-sig0),'--','tag','curves','DisplayName','\mu\pm\sigma')
handles.ax1.ColorOrderIndex=1; plot(handles.ax1,cav0,exp(lny0+sig0),'--','tag','curves','HandleVisibility','off')
Leg=legend(handles.ax1);
Leg.Box      = 'off';
Leg.Location = 'NorthWest';

c2 = uicontextmenu;
uimenu(c2,'Label','Copy data','Callback'  ,{@data2clipboard_uimenu,num2cell([cav0(:),exp([lny0-sig0;lny0;lny0+sig0])'])});
set(handles.ax1,'uicontextmenu',c2);
% 
delete(findall(handles.ax2,'tag','curves'))
handles.ax2.ColorOrderIndex=1;
plot(handles.ax2,s,Y,'tag','curves','HandleVisibility','off')

c3 = uicontextmenu;
uimenu(c3,'Label','Copy data','Callback'  ,{@data2clipboard_uimenu,num2cell([s(:),Y(:)])});
set(handles.ax2,'uicontextmenu',c3);

% ---------------------- building properties ------------------------------

function B_Callback(hObject, eventdata, handles) %#ok<*DEFNU>
plotSettle(handles)

function L_Callback(hObject, eventdata, handles)
plotSettle(handles)

function Df_Callback(hObject, eventdata, handles)
plotSettle(handles)

function Q_Callback(hObject, eventdata, handles)
if isempty(hObject.String)
    handles.numStories.Enable='on';
else
    handles.numStories.Enable='off';
end
plotSettle(handles)

function bType_Callback(hObject, eventdata, handles)
plotSettle(handles)

function wt_Callback(hObject, eventdata, handles)
plotSettle(handles)
% ------------------------ soil deposit parameters ------------------------

function LPC_Callback(hObject, eventdata, handles)
plotSettle(handles)

function th1_Callback(hObject, eventdata, handles)
plotSettle(handles)

function th2_Callback(hObject, eventdata, handles)
plotSettle(handles)

function N1_Callback(hObject, eventdata, handles)
plotSettle(handles)

function meth_Callback(hObject, eventdata, handles)

switch hObject.Value
    case 1, handles.table1.ColumnName{1}='N1_60';
    case 2, handles.table1.ColumnName{1}='qc1N';
end
plotSettle(handles)

function N2_Callback(hObject, eventdata, handles)
plotSettle(handles)

function LBS_Callback(hObject, eventdata, handles) %#ok<*INUSD>

function HL_Callback(hObject, eventdata, handles)

% --------------- Settlement / Tilt Model ---------------------------------
function STpop_Callback(hObject, eventdata, handles)
handles=popSettleModel(handles);
plotSettle(handles)
guidata(hObject,handles)

function e1_Callback(hObject, eventdata, handles)
plotSettle(handles)

function e2_Callback(hObject, eventdata, handles)
plotSettle(handles)

function e3_Callback(hObject, eventdata, handles)
plotSettle(handles)

% --------------- Other ---------------------------------------------------
function DiscretizeMButtonbutton_Callback(hObject, eventdata, handles)
switch handles.plotmodeax2
    case 'pdf'  , handles.plotmodeax2='cdf';
    case 'cdf'  , handles.plotmodeax2='ccdf';
    case 'ccdf' , handles.plotmodeax2='pdf';
end
plotSettle(handles)
guidata(hObject,handles)

function undock1_Callback(hObject, eventdata, handles)
figure2clipboard_uimenu(hObject, eventdata,handles.ax1)

function undock2_Callback(hObject, eventdata, handles)
figure2clipboard_uimenu(hObject, eventdata,handles.ax2)

function Book_Callback(hObject, eventdata, handles)
val     = handles.STpop.Value;
methods = pshatoolbox_methods(6,val);
if ~isempty(methods.ref)
    try
        web(methods.ref,'-browser')
    catch
    end
end



function e4_Callback(hObject, eventdata, handles)
