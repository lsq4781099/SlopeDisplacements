function varargout = MR_explorer(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @MR_explorer_OpeningFcn, ...
    'gui_OutputFcn',  @MR_explorer_OutputFcn, ...
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

function MR_explorer_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
handles.Exit_button.CData=double(imread('Exit.jpg'))/255;
handles.openbook.CData=double(imread('book_open.jpg'))/255;
handles.AxisScale.CData=double(imresize(imread('Ruler.jpg'),[20 20]))/255;
handles.gridmanager.CData=double(imread('Grid.jpg'))/255;

handles.xlabel=xlabel(handles.ax1,'Magnitude','fontsize',10);
handles.ylabel=ylabel(handles.ax1,'Magnitude pdf','fontsize',10);

methods = pshatoolbox_methods(2);
handles.MRselect.String = {methods.label};

handles.ax1.Box='on';
handles.ax1.Color=[1 1 1];
handles.ax1.XGrid='on';
handles.ax1.YGrid='on';
handles.ax1.NextPlot='add';

% initialize data to fill plot
handles.paramlist = cell(0,2);
handles.ptrs      = zeros(0,15);
ch  = get(handles.panel2,'children'); tag  = char(ch.Tag);
handles.text = flipud(find(ismember(tag(:,1),'t')));
handles.edit = flipud(find(ismember(tag(:,1),'e')));
handles.selectedrow = [];

set(ch(handles.text),'Visible','off')
set(ch(handles.edit),'Visible','off');
handles.uitable1.Data=cell(0,2);

if nargin==10
    handles.mrr1 = varargin{1};
    handles.mrr2 = varargin{2};
    handles.mrr3 = varargin{3};
    handles.mrr4 = varargin{4};
    handles.mrr5 = varargin{5};
    Nmrr  = varargin{6};
    handles.mu  = varargin{7}.ShearModulus;
    Nmscl = length(handles.mrr1);
    mscl.source= cell(0,1);
    mscl.ptr   = zeros(0,3);
    fhand = cell(0,1);
    for i=1:Nmscl
        mscl.source = [mscl.source;handles.mrr1(i).source;handles.mrr2(i).source;handles.mrr3(i).source;handles.mrr4(i).source;handles.mrr5(i).source];
        n1 = Nmrr(1,i);
        n2 = Nmrr(2,i);
        n3 = Nmrr(3,i);
        n4 = Nmrr(4,i);
        n5 = Nmrr(5,i);
        mscl.ptr = [mscl.ptr;
            ones(n1,1)*i,ones(n1,1)*1,(1:n1)';
            ones(n2,1)*i,ones(n2,1)*2,(1:n2)';
            ones(n3,1)*i,ones(n3,1)*3,(1:n3)';
            ones(n4,1)*i,ones(n4,1)*4,(1:n4)';
            ones(n5,1)*i,ones(n5,1)*5,(1:n5)'];
        
        fhand=[fhand;
            repmat({'delta'}     ,n1,1);
            repmat({'truncexp'}  ,n2,1);
            repmat({'truncnorm'} ,n3,1);
            repmat({'yc1985'}    ,n4,1);
            repmat({'magtable'}  ,n5,1);
            ]; %#ok<AGROW>
    end
    handles.mscl=mscl;
    handles.uitable1.Data=[mscl.source,fhand];
    handles = CellSelectAction(handles,1);
    handles.AddNew.Enable='off';
else
    handles = mMRdefault(handles,ch(handles.text),ch(handles.edit));
    plotmrmodel(handles);
end

axis(handles.ax1,'auto')
guidata(hObject, handles);
% uiwait(handles.figure1);

function varargout = MR_explorer_OutputFcn(hObject, eventdata, handles)  %#ok<*INUSL>
varargout{1} = handles.output;

% ------------ key funtions ---------------------------------------------------
function MRselect_Callback(hObject, eventdata, handles) %#ok<*INUSD,*DEFNU>

ch=get(handles.panel2,'children');
set(ch(handles.text),'Visible','off')
set(ch(handles.edit),'Visible','off','Style','edit');

handles = mMRdefault(handles,ch(handles.text),ch(handles.edit));
plotmrmodel(handles);
guidata(hObject,handles)

function Exit_button_Callback(hObject, eventdata, handles)
close(handles.figure1)

function AxisScale_Callback(hObject, eventdata, handles)
XYSCALE=[handles.ax1.XScale,handles.ax1.YScale];
switch XYSCALE
    case 'linearlinear', handles.ax1.XScale='log';
    case 'loglinear',    handles.ax1.XScale='linear'; handles.ax1.YScale='log';
    case 'linearlog',    handles.ax1.XScale='log';
    case 'loglog',       handles.ax1.XScale='linear'; handles.ax1.YScale='linear';
end

function AutoScale_Callback(hObject, eventdata, handles)
switch hObject.Value
    case 0
        prompt={'Xscale','Yscale','XLim','YLim'};
        name='Manual Scale'; numlines=1;
        XL = handles.ax1.XLim; strX = [num2str(XL(1)),' ',num2str(XL(2))];
        YL = handles.ax1.YLim; strY = [num2str(YL(1)),' ',num2str(YL(2))];
        defaultanswer={handles.ax1.XScale,handles.ax1.YScale,strX,strY};
        answer=inputdlg(prompt,name,numlines,defaultanswer);
        if isempty(answer)
            return
        end
        XLIM = regexp(answer{3},'\s','split'); YLIM = regexp(answer{4},'\s','split');
        handles.ax1.XScale=answer{1};
        handles.ax1.YScale=answer{2};
        handles.ax1.XLim=[eval(XLIM{1}) eval(XLIM{2})];
        handles.ax1.YLim=[eval(YLIM{1}) eval(YLIM{2})];
    case 1
        axis(handles.ax1,'auto')
end
guidata(hObject,handles)

function HoldPlot_Callback(hObject, eventdata, handles)

function File_Callback(hObject, eventdata, handles)

function Library_Callback(hObject, eventdata, handles)

function AddNew_Callback(hObject, eventdata, handles)
guidata(hObject,handles)

function UpdateModel_Callback(hObject, eventdata, handles)
guidata(hObject,handles)

function LoadList_Callback(hObject, eventdata, handles)

function RemoveSelection_Callback(hObject, eventdata, handles)
guidata(hObject,handles)

% ------------  edit boxes
function e1_Callback(hObject, eventdata, handles)
plotmrmodel(handles);
guidata(hObject,handles)

function e2_Callback(hObject, eventdata, handles)

plotmrmodel(handles);
guidata(hObject,handles)

function e3_Callback(hObject, eventdata, handles)
plotmrmodel(handles);
guidata(hObject,handles)

function e4_Callback(hObject, eventdata, handles)
plotmrmodel(handles);
guidata(hObject,handles)

function e5_Callback(hObject, eventdata, handles)
plotmrmodel(handles);
guidata(hObject,handles)

function e6_Callback(hObject, eventdata, handles)

plotmrmodel(handles);
guidata(hObject,handles)

function e7_Callback(hObject, eventdata, handles)
plotmrmodel(handles);
guidata(hObject,handles)

function e8_Callback(hObject, eventdata, handles)
plotmrmodel(handles);
guidata(hObject,handles)

function e9_Callback(hObject, eventdata, handles)

switch handles.MRselect.String{handles.MRselect.Value}
    case 'Abrahamson Silva 2008 - NGA' % set default Ztop
        Vs30 = str2double(handles.e9.String);
        handles.e8.String=sprintf('%4.4g',Z10_default_AS08_NGA(Vs30));
end

plotmrmodel(handles);
guidata(hObject,handles)

function e10_Callback(hObject, eventdata, handles)
switch handles.MRselect.String{handles.MRselect.Value}
    case 'Chiou Youngs 2008 - NGA' % set default Z1.0
        Vs30 = str2double(handles.e10.String);
        %handles.e7.String=sprintf('%4.4g',exp(28.5-3.82/8*log(Vs30^8+378.8^8)));
        handles.e7.String=sprintf('%4.4g',Z10_default_AS08_NGA(Vs30));
end

plotmrmodel(handles);
guidata(hObject,handles)

function e11_Callback(hObject, eventdata, handles)
plotmrmodel(handles);
guidata(hObject,handles)

function e12_Callback(hObject, eventdata, handles)
plotmrmodel(handles);
guidata(hObject,handles)

function e13_Callback(hObject, eventdata, handles)
plotmrmodel(handles);
guidata(hObject,handles)

function e14_Callback(hObject, eventdata, handles)
plotmrmodel(handles);
guidata(hObject,handles)

function gridmanager_Callback(hObject, eventdata, handles)
switch [handles.ax1.XGrid,handles.ax1.XMinorGrid]
    case 'offoff'
        handles.ax1.XGrid     ='on';  handles.ax1.XMinorGrid='on';
        handles.ax1.YGrid     ='on';  handles.ax1.YMinorGrid='on';
    case 'onon'
        handles.ax1.XGrid     ='on';  handles.ax1.XMinorGrid='off';
        handles.ax1.YGrid     ='on';  handles.ax1.YMinorGrid='off';
    case 'onoff'
        handles.ax1.XGrid     ='off';  handles.ax1.XMinorGrid='on';
        handles.ax1.YGrid     ='off';  handles.ax1.YMinorGrid='on';
    case 'offon'
        handles.ax1.XGrid     ='off';  handles.ax1.XMinorGrid='off';
        handles.ax1.YGrid     ='off';  handles.ax1.YMinorGrid='off';
end

% --------------------------------------------------------------------
function grouplist_Callback(hObject, eventdata, handles)

function grouplist_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function uitable1_CellSelectionCallback(hObject, eventdata, handles)
if isempty(eventdata.Indices)
    return
end
ind     = eventdata.Indices(1);
handles = CellSelectAction(handles,ind);
guidata(hObject,handles)

function[handles]=CellSelectAction(handles,ind)

ms     = handles.mscl.ptr(ind,1);
mstype = handles.mscl.ptr(ind,2);
mspos  = handles.mscl.ptr(ind,3);

handles.MRselect.Value = mstype;
ch=get(handles.panel2,'children');
set(ch(handles.text),'Visible','off')
set(ch(handles.edit),'Visible','off','Style','edit');
handles = mMRdefault(handles,ch(handles.text),ch(handles.edit));

switch mstype
    case 1
        num=handles.mrr1(ms).num(mspos,:);
        [~,~,meanMo]=delta(5:8,num(3));
        if num(1)==2
            A = 100*(100000^2); %100 km2
            S = num(2);
            NMmin = handles.mu*A*S/meanMo;
        else
            NMmin = num(2);
        end
        
        handles.e1.String=sprintf('%g',NMmin);
        handles.e2.String=sprintf('%g',num(3));
        if num(1)==2
            handles.e3.String=sprintf('%g',S);
            handles.e4.String=sprintf('%g',meanMo);
        else
            handles.e3.String='-';
            handles.e4.String='-';
        end
    case 2
        num=handles.mrr2(ms).num(mspos,:);
        [~,~,meanMo]=truncexp(5:8,num(3:5));
        if num(1)==2
            A = 100*(100000^2); %100 km2
            S = num(2);
            NMmin = handles.mu*A*S/meanMo;
        else
            NMmin = num(2);
        end
        
        handles.e1.String=sprintf('%g',NMmin);
        handles.e2.String=sprintf('%g',num(3));
        handles.e3.String=sprintf('%g',num(4));
        handles.e4.String=sprintf('%g',num(5));
        if num(1)==2
            handles.e5.String=sprintf('%g',S);
            handles.e6.String=sprintf('%g',meanMo);
        else
            handles.e5.String='-';
            handles.e6.String='-';
        end
        
    case 3
        num=handles.mrr3(ms).num(mspos,:);
        [~,~,meanMo]=truncnorm(5:8,num(3:6));
        if num(1)==2
            A = 100*(100000^2); %100 km2
            S = num(2);
            NMmin = handles.mu*A*S/meanMo;
        else
            NMmin = num(2);
        end
        handles.e1.String=sprintf('%g',NMmin);
        handles.e2.String=sprintf('%g',num(3));
        handles.e3.String=sprintf('%g',num(4));
        handles.e4.String=sprintf('%g',num(5));
        handles.e5.String=sprintf('%g',num(6));
        if num(1)==2
            handles.e6.String=sprintf('%g',S);
            handles.e7.String=sprintf('%g',meanMo);
        else
            handles.e6.String='-';
            handles.e7.String='-';
        end
        
        
    case 4
        num=handles.mrr4(ms).num(mspos,:);
        [~,~,meanMo]=yc1985(5:8,num(3:5));
        if num(1)==2
            A = 100*(100000^2); %100 km2
            S = num(2);
            NMmin = handles.mu*A*S/meanMo;
        else
            NMmin = num(2);
        end
        handles.e1.String=sprintf('%g',NMmin);
        handles.e2.String=sprintf('%g',num(3));
        handles.e3.String=sprintf('%g',num(4));
        handles.e4.String=sprintf('%g',num(5));
        if num(1)==2
            handles.e5.String=sprintf('%g',S);
            handles.e6.String=sprintf('%g',meanMo);
        else
            handles.e5.String='-';
            handles.e6.String='-';
        end
    case 5
        num=handles.mrr5(ms).num{mspos};
        NMmin = num(2);
        Mmin  = num(3);
        dM    = num(4);
        occr  = num(5:end);
        handles.e1.String=sprintf('%g',Mmin);
        handles.e2.String=sprintf('%g',dM);
        handles.e3.String=sprintf(repmat('%g ',1,length(occr)),occr');
end

plotmrmodel(handles);
handles.selectedrow=ind;

function uitable1_CellEditCallback(hObject, eventdata, handles)

guidata(hObject,handles)

% --------------------------------------------------------------------
function UndockAxis_Callback(hObject, eventdata, handles)
figure2clipboard_uimenu(hObject, eventdata,handles.ax1)

function uibuttongroup1_SelectionChangedFcn(hObject, eventdata, handles)
plotmrmodel(handles);
