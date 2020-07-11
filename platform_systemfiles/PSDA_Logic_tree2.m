function varargout = PSDA_Logic_tree2(varargin)
%#ok<*AGROW>
%#ok<*INUSD>
%#ok<*INUSL>
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @PSDA_Logic_tree2_OpeningFcn, ...
    'gui_OutputFcn',  @PSDA_Logic_tree2_OutputFcn, ...
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

function PSDA_Logic_tree2_OpeningFcn(hObject, ~, handles, varargin)
handles.output = [];
handles.sys       = varargin{1};
h                 = varargin{2}; 
T1                = varargin{3};
T2                = varargin{4};
T3                = varargin{5};
AnalysisType      = varargin{6};
handles.paramPSDA = varargin{7};

handles.sitepop.Value  = 1; 
handles.sitepop.String = h.id;

switch AnalysisType
    case 'PBPA'  , handles.radiobutton1.Value=1;
    case 'FPBPA' , handles.radiobutton2.Value=1;
end

label    = {handles.sys.SMLIB.id};
isGMMReg = horzcat(handles.sys.SMLIB.isregular);

handles.table3.ColumnFormat{2}= label(isGMMReg);
handles.table3.ColumnFormat{3}= label(isGMMReg);
handles.table3.ColumnFormat{4}= label(isGMMReg);

% seismic hazard
handles.table1.Data=T1;

% slope parameters
[~,kyptr]=intersect(h.param,{'ky','covky'},'stable');
[~,Tsptr]=intersect(h.param,{'Ts','covTs'},'stable');
handles.ky = h.value(:,kyptr);
handles.Ts = h.value(:,Tsptr);

val = 1;
handles.Ts_mean.String = sprintf('%g',handles.Ts(val,1));
handles.Ts_cov.String  = sprintf('%g',handles.Ts(val,2));
handles.Ts_nsta.String = sprintf('%g',handles.paramPSDA.Tssamples);
handles.ky_mean.String = sprintf('%g',handles.ky(val,1));
handles.ky_cov.String  = sprintf('%g',handles.ky(val,2));
handles.ky_nsta.String = sprintf('%g',handles.paramPSDA.kysamples);

handles.table2.Data    = T2;

% displacement model
handles.table3.Data = T3;
handles.current1    = zeros(0,2);
handles.current2    = zeros(0,2);
handles.current3    = zeros(0,2);

% axis properties
handles.ax1.Box='off';
handles.ax1.XTick=[];
handles.ax1.YTick=[];
handles.ax1.XLim=[0 9];
handles.ax1.YLim=[-1.1 1.3];
handles.ax1.NextPlot='add';
text(handles.ax1,2,1.3,'Seismic Hazard','horizontalalignment','right','fontsize',8)
text(handles.ax1,5,1.3,'Slope Parameters','horizontalalignment','right','fontsize',8)
text(handles.ax1,8,1.3,'Displacement Model','horizontalalignment','right','fontsize',8)

switch handles.radiobutton1.Value
    case 1, plotPSDALogicTree(handles);
    case 0, plotPSDALogicTree(handles);
end
axis(handles.ax1,'off')
guidata(hObject, handles);
uiwait(handles.figure1);

function varargout = PSDA_Logic_tree2_OutputFcn(hObject, eventdata, handles)

% normalizacion de los pesos, just in case...
% -------------------------------------------------------------------------
w1 = cell2mat(handles.table1.Data(:,2)); w1 = w1/sum(w1); handles.table1.Data(:,2)=num2cell(w1);
w2 = cell2mat(handles.table2.Data(:,4)); w2 = w2/sum(w2); handles.table2.Data(:,4)=num2cell(w2);
w3 = cell2mat(handles.table3.Data(:,5)); w3 = w3/sum(w3); handles.table3.Data(:,5)=num2cell(w3);
% -------------------------------------------------------------------------

switch handles.radiobutton1.Value
    case 1
        varargout{1}={'Branch0',1};
    case 0
        varargout{1}=handles.table1.Data;
end
varargout{2}=handles.table2.Data;
varargout{3}=handles.table3.Data;
varargout{4}=handles.uibuttongroup2.SelectedObject.String;
varargout{5}=handles.paramPSDA;
delete(handles.figure1)

function plotPSDALogicTree(handles)
delete(findall(handles.ax1,'tag','logic'));
handles.ax1.ColorOrderIndex=1;
N1 = size(handles.table1.Data,1);
N2 = size(handles.table2.Data,1);
N3 = size(handles.table3.Data,1);

x1  = [0 1 2];  LT1 = buildLTLine(x1,N1);
x2  = x1+3;     LT2 = buildLTLine(x2,N2);
x3  = x1+6;     LT3 = buildLTLine(x3,N3);

plot(handles.ax1,LT1(:,1),LT1(:,2),'.-','tag','logic'),handles.ax1.ColorOrderIndex=1;
plot(handles.ax1,LT2(:,1),LT2(:,2),'.-','tag','logic'),handles.ax1.ColorOrderIndex=1;
plot(handles.ax1,LT3(:,1),LT3(:,2),'.-','tag','logic'),handles.ax1.ColorOrderIndex=1;

% c    = uicontextmenu;
% uimenu(c,'Label','Undock','Callback',           {@figure2clipboard_uimenu,handles.ax1,''});
% set(handles.ax1,'uicontextmenu',c);

function LT = buildLTLine(x,N)
if N==1
    y  = 0;
else
    y  = linspace(-1,1,N);
end
LT = [x(1) 0;NaN NaN];
for i=1:N
    nLT = [x(1) 0;x(2) y(i);x(3) y(i);NaN NaN];
    LT  = [LT;nLT]; 
end

function Normalize2_Callback(hObject, eventdata, handles) %#ok<*DEFNU>

w = handles.table2.Data(:,4);
w = cell2mat(w);
handles.table2.Data(:,4)=num2cell(w/sum(w));
guidata(hObject,handles)

function Ts_mean_Callback(hObject, eventdata, handles)

function Ts_cov_Callback(hObject, eventdata, handles)

function Ts_nsta_Callback(hObject, eventdata, handles)
handles.paramPSDA.Tssamples = max(str2double(hObject.String),0);
val = handles.sitepop.Value;
handles.table2.Data = buildPSDA_T2(handles.paramPSDA,handles.ky(val,:),handles.Ts(val,:));
plotPSDALogicTree(handles)
guidata(hObject,handles)

function ky_mean_Callback(hObject, eventdata, handles)

function ky_cov_Callback(hObject, eventdata, handles)

function ky_nsta_Callback(hObject, eventdata, handles)
handles.paramPSDA.kysamples = max(str2double(hObject.String),0);
val = handles.sitepop.Value;
handles.table2.Data = buildPSDA_T2(handles.paramPSDA,handles.ky(val,:),handles.Ts(val,:));
plotPSDALogicTree(handles)
guidata(hObject,handles)

function Normalize3_Callback(hObject, eventdata, handles)

w = handles.table3.Data(:,5);
wsum = sum(cell2mat(w));
for i=1:length(w)
    if ~isempty(w{i})
        w{i}=w{i}/wsum;
    end
end
handles.table3.Data(1:length(w),5)=w;
guidata(hObject,handles)

function New_Disp_Callback(hObject, eventdata, handles)
handles.table3.Data(end+1,:)=handles.table3.Data(end,:);
handles.current3=size(handles.table3.Data);

for i=1:size(handles.table3.Data,1)
    handles.table3.Data{i,1}=sprintf('slope%g',i);
end
plotPSDALogicTree(handles)
guidata(hObject,handles)

function Delete_Disp_Callback(hObject, eventdata, handles)

if ~isempty(handles.current3)
    if size(handles.table3.Data,1)>1
        ind = handles.current3(1);
        handles.table3.Data(ind,:)=[];
        plotPSDALogicTree(handles)
        handles.current3=zeros(0,2);
    end
end

for i=1:size(handles.table3.Data,1)
    handles.table3.Data{i,1}=sprintf('disp%g',i);
end

guidata(hObject,handles)

function table1_CellEditCallback(hObject, eventdata, handles)
handles.current1=eventdata.Indices;
guidata(hObject,handles);

function figure1_CloseRequestFcn(hObject, eventdata, handles) 
if isequal(get(hObject,'waitstatus'),'waiting')
    uiresume(hObject);
else
    delete(hObject);
end
% delete(hObject);

function uibuttongroup2_SelectionChangedFcn(hObject, eventdata, handles)

switch hObject.String
    case 'PBPA'  
        handles.table1.Data={'Branch0',1};
        plotPSDALogicTree(handles)
    case 'FPBPA' 
        W = handles.sys.branch(:,4);
        Nb = length(W);
        handles.table1.Data=[compose('Branch%i',(1:Nb)'),num2cell(W(:))];
        plotPSDALogicTree(handles)
end

function sitepop_Callback(hObject, eventdata, handles)

val = hObject.Value;
handles.Ts_mean.String = sprintf('%g',handles.Ts(val,1));
handles.Ts_cov.String  = sprintf('%g',handles.Ts(val,2));
handles.Ts_nsta.String = sprintf('%g',handles.paramPSDA.Tssamples);
handles.ky_mean.String = sprintf('%g',handles.ky(val,1));
handles.ky_cov.String  = sprintf('%g',handles.ky(val,2));
handles.ky_nsta.String = sprintf('%g',handles.paramPSDA.kysamples);

handles.table2.Data = buildPSDA_T2(handles.paramPSDA,handles.ky(val,:),handles.Ts(val,:));  
guidata(hObject,handles)
