function varargout = LIBSLogicTree(varargin)
%#ok<*INUSD>
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @LIBSLogicTree_OpeningFcn, ...
    'gui_OutputFcn',  @LIBSLogicTree_OutputFcn, ...
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

function LIBSLogicTree_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>
handles.output = [];
SETLIB     = varargin{1};
T1         = varargin{2};
T2         = varargin{3};
T3         = varargin{4};
optlib     = varargin{5};
weightHAZ  = varargin{6};

handles.table3.ColumnFormat{2}={SETLIB.label};
handles.table3.ColumnFormat{3}={SETLIB.label};
handles.table3.ColumnFormat{4}={SETLIB.label};

switch optlib.analysis
    case 'PBPA'  , handles.radiobutton1.Value=1;
    case 'FPBPA' , handles.radiobutton2.Value=1;
end

handles.table1.Data = T1;
handles.table2.Data = T2(:,2:5);
handles.table3.Data = T3;
handles.nQ.String   = optlib.nQ;
handles.nPGA.String = optlib.nPGA;
handles.wHL1.String = optlib.wHL(1);
handles.wHL2.String = optlib.wHL(2);
handles.RetPeriod.String = optlib.RetPeriod;
handles.optlib      = optlib;
handles.weightsHAZ  = weightHAZ;

plotLIBSLogicTree(handles)

% axis properties
axis(handles.ax1,'off');
handles.ax1.Box='off';
handles.ax1.XTick=[];
handles.ax1.YTick=[];
text(handles.ax1,2,1.3,'Seismic Hazard Model','horizontalalignment','right','fontsize',8)
text(handles.ax1,5,1.3,'Site and Building Model','horizontalalignment','right','fontsize',8)
text(handles.ax1,8,1.3,'Settlement Model','horizontalalignment','right','fontsize',8)

handles.current1    = zeros(0,2);
handles.current2    = zeros(0,2);
handles.current3    = zeros(0,2);

guidata(hObject, handles);
uiwait(handles.figure1);

function varargout = LIBSLogicTree_OutputFcn(hObject, eventdata, handles)

nT2             = size(handles.table2.Data,1);
optlib          = handles.optlib;
optlib.analysis = handles.uibuttongroup1.SelectedObject.String;
optlib.nQ       = str2double(handles.nQ.String);
optlib.nPGA     = str2double(handles.nPGA.String);
optlib.wHL      = [str2double(handles.wHL1.String),str2double(handles.wHL2.String)];
optlib.RetPeriod = str2double(handles.RetPeriod.String);
varargout{1} = handles.table1.Data;
varargout{2} = [compose('SP%i',1:nT2)',handles.table2.Data];
varargout{3} = handles.table3.Data;
varargout{4} = optlib;
delete(handles.figure1)

function plotLIBSLogicTree(handles)
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

function LT = buildLTLine(x,N)
if N==1
    y  = 0;
else
    y  = linspace(-1,1,N);
end
LT = [x(1) 0;NaN NaN];
for i=1:N
    nLT = [x(1) 0;x(2) y(i);x(3) y(i);NaN NaN];
    LT  = [LT;nLT]; %#ok<*AGROW>
end

function table1_CellSelectionCallback(hObject, eventdata, handles) %#ok<*DEFNU>

if ~isempty(eventdata.Indices)
    ind = eventdata.Indices(1);
    handles.model_t1.String=handles.model(ind).id1;
    handles.model_t2.String=handles.model(ind).id2;
    handles.model_t3.String=handles.model(ind).id3;
end
guidata(hObject,handles)

function Normalize1_Callback(hObject, eventdata, handles)

w = handles.table1.Data(:,2);
w = cell2mat(w);
handles.table1.Data(:,2)=num2cell(w/sum(w));
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
plotLIBSLogicTree(handles)
guidata(hObject,handles)

function Delete_Disp_Callback(hObject, eventdata, handles)

if ~isempty(handles.current3)
    if size(handles.table3.Data,1)>1
        ind = handles.current3(1);
        handles.table3.Data(ind,:)=[];
        plotLIBSLogicTree(handles)
        handles.current3=zeros(0,2);
    end
end
for i=1:size(handles.table3.Data,1)
    handles.table3.Data{i,1}=sprintf('ST%i',i);
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

function table3_CellEditCallback(hObject, eventdata, handles)

function nQ_Callback(hObject, eventdata, handles)
nQ   = max(str2double(handles.nQ.String),1);
nPGA = str2double(handles.nPGA.String);
wHL  = [str2double(handles.wHL1.String),str2double(handles.wHL2.String)];
T2   = T2settle(nQ,nPGA,wHL);
handles.table2.Data=num2cell(T2);
plotLIBSLogicTree(handles)

function nPGA_Callback(hObject, eventdata, handles)
nQ   = str2double(handles.nQ.String);
nPGA = max(str2double(handles.nPGA.String),1);
wHL  = [str2double(handles.wHL1.String),str2double(handles.wHL2.String)];
T2   = T2settle(nQ,nPGA,wHL);
handles.table2.Data=num2cell(T2);
plotLIBSLogicTree(handles)

function wHL1_Callback(hObject, eventdata, handles)

wHL1 = min(max(str2double(hObject.String),0),1);
handles.wHL1.String=num2str(wHL1);
handles.wHL2.String=num2str(1-wHL1);
nQ   = str2double(handles.nQ.String);
nPGA = str2double(handles.nPGA.String);
wHL  = [str2double(handles.wHL1.String),str2double(handles.wHL2.String)];
T2   = T2settle(nQ,nPGA,wHL);
handles.table2.Data=num2cell(T2);
plotLIBSLogicTree(handles)

function wHL2_Callback(hObject, eventdata, handles)
wHL2 = min(max(str2double(hObject.String),0),1);
handles.wHL1.String=num2str(1-wHL2);
handles.wHL2.String=num2str(wHL2);
nQ   = str2double(handles.nQ.String);
nPGA = str2double(handles.nPGA.String);
wHL  = [str2double(handles.wHL1.String),str2double(handles.wHL2.String)];
T2   = T2settle(nQ,nPGA,wHL);
handles.table2.Data=num2cell(T2);
plotLIBSLogicTree(handles)

function T2=T2settle(nQ,nPGA,wHL)

wHL=wHL(:);
% wHL(wHL==0)=[];
[~,~,wQ]   = trlognpdf_psda([1 0.2 nQ]);
[~,~,wLBS] = trlognpdf_psda([1 0.2 nPGA]);
[II,JJ,KK]=meshgrid(1:nQ,1:nPGA,1:length(wHL)); II=II(:); JJ=JJ(:);KK=KK(:);
T2=[wQ(II),wLBS(JJ),wHL(KK),wQ(II).*wLBS(JJ).*wHL(KK)];

function uibuttongroup1_SelectionChangedFcn(hObject, eventdata, handles)

switch hObject.String
    case 'PBPA'  
        handles.table1.Data={'Haz0',1};
        plotLIBSLogicTree(handles)
    case 'FPBPA' 
        W = handles.weightsHAZ;
        Nb = length(W);
        handles.table1.Data=[compose('Haz%i',(1:Nb)'),num2cell(W(:))];
        plotLIBSLogicTree(handles)
end

function table3_CellSelectionCallback(hObject, eventdata, handles)
handles.current3=eventdata.Indices;
guidata(hObject,handles)

function RetPeriod_Callback(hObject, eventdata, handles)
