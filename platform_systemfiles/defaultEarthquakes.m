function varargout = defaultEarthquakes(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @defaultEarthquakes_OpeningFcn, ...
                   'gui_OutputFcn',  @defaultEarthquakes_OutputFcn, ...
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

function defaultEarthquakes_OpeningFcn(hObject, eventdata, handles, varargin) 

if nargin==3
    SUB  = createObj('interfaceEQ');
    SC   = createObj('crustalEQ');
    SITE = createObj('siteGMM');
else
    SUB  = varargin{1};
    SC   = varargin{2};
    SITE = varargin{3};
    
end

handles.popupmenu1.Value = SUB.Mech;
handles.edit1.String  = SUB.Mag;
handles.edit2.String  = SUB.Ztor;
handles.edit3.String  = SUB.Rrup;
handles.edit4.String  = SUB.Rx;
handles.edit5.String  = SUB.Rhyp;
handles.edit6.String  = SUB.Zhyp;

handles.popupmenu2.Value=SC.Mech;
handles.edit11.String = SC.Mag; 
handles.edit12.String = SC.dip; 
handles.edit13.String = SC.W;   
handles.edit14.String = SC.Zbot;
handles.edit15.String = SC.Ztor;
handles.edit16.String = SC.Rx;  
handles.edit17.String = SC.Rrup;
handles.edit18.String = SC.Rjb; 
handles.edit19.String = SC.Rhyp;
handles.edit20.String = SC.Zhyp;
handles.edit21.String = SC.Ry0; 
handles.popupmenu3.Value=SC.region;

handles.VS30.String   = SITE.VS30;
handles.f0.String     = SITE.f0;
handles.Z10.String    = SITE.Z10;
handles.Z25.String    = SITE.Z25;

switch SITE.Idini
    case 1, handles.radiobutton5.Value  = 1;
    case 2, handles.radiobutton6.Value  = 1;
    case 3, handles.radiobutton7.Value  = 1;
    case 4, handles.radiobutton8.Value  = 1;
    case 5, handles.radiobutton9.Value  = 1;
    case 6, handles.radiobutton10.Value = 1;
end

switch SITE.SGS
    case 1, handles.radiobutton2.Value  = 1;
    case 2, handles.radiobutton3.Value  = 1;
    case 3, handles.radiobutton4.Value  = 1;
end

% stores default data
handles.defaultSUB    = SUB;
handles.defaultSC     = SC;
handles.defaultSITE   = SITE;

guidata(hObject, handles);
uiwait(handles.figure1);

function varargout = defaultEarthquakes_OutputFcn(hObject, eventdata, handles)  %#ok<*INUSL>

% subduiction default parameters
SUB.Mech = handles.popupmenu1.Value;
SUB.Mag  = str2double(handles.edit1.String) ;
SUB.Ztor = str2double(handles.edit2.String) ;
SUB.Rrup = str2double(handles.edit3.String) ;
SUB.Rx   = str2double(handles.edit4.String) ;
SUB.Rhyp = str2double(handles.edit5.String) ;
SUB.Zhyp = str2double(handles.edit6.String) ;

Nsamples    = 40;
SUB.rpMag   = SUB.Mag*ones(Nsamples,1);
SUB.rpZtor  = SUB.Ztor*ones(Nsamples,1);
SUB.rpRrup  = logsp(SUB.rpZtor(1),400,Nsamples)';
SUB.rpRx    = sqrt(max(SUB.rpRrup.^2-SUB.rpZtor.^2,0));
SUB.rpRhyp  = SUB.rpRrup;
SUB.rpZhyp  = SUB.rpZtor;


% shallow crustal default parameters
SC.Mech  = handles.popupmenu2.Value;
SC.Mag   = str2double(handles.edit11.String);
SC.dip   = str2double(handles.edit12.String);
SC.W     = str2double(handles.edit13.String);
SC.Zbot  = str2double(handles.edit14.String);
SC.Ztor  = str2double(handles.edit15.String);
SC.Rx    = str2double(handles.edit16.String);
SC.Rrup  = str2double(handles.edit17.String);
SC.Rjb   = str2double(handles.edit18.String);
SC.Rhyp  = str2double(handles.edit19.String);
SC.Zhyp  = str2double(handles.edit20.String);
SC.Ry0   = str2double(handles.edit21.String);
SC.Ry0   = str2double(handles.edit21.String);
SC.region= handles.popupmenu3.Value;


Nsamples  = 40;
SC2.Mag   = SC.Mag*ones(Nsamples,1);
SC2.dip   = SC.dip;
SC2.W     = SC.W;
SC2.Zbot  = SC.Zbot;
SC2.Ztor  = SC.Ztor*ones(Nsamples,1);
SC2.Rx    = [0;logsp(1,400,Nsamples-1)'];
SC2       = getRrupRjb(SC2);
SC2.Ry0   = 0*ones(Nsamples,1);
SC.rpMag  = SC2.Mag;
SC.rpZtor = SC2.Ztor;
SC.rpRx   = SC2.Rx;
SC.rpRrup = SC2.Rrup;
SC.rpRjb  = SC2.Rjb;
SC.rpRhyp = SC2.Rhyp;
SC.rpZhyp = SC2.Zhyp;
SC.rpRy0  = SC2.Ry0;

SITE.VS30  = str2double(handles.VS30.String);
SITE.f0    = str2double(handles.f0.String);
SITE.Z10   = str2double(handles.Z10.String);
SITE.Z25   = str2double(handles.Z25.String);

switch handles.uibuttongroup2.SelectedObject.String
    case 'sI'  , SITE.Idini = 1;
    case 'sII' , SITE.Idini = 2;
    case 'sIII', SITE.Idini = 3;
    case 'sIV' , SITE.Idini = 4;
    case 'sV'  , SITE.Idini = 5;
    case 'sVI' , SITE.Idini = 6;
end

switch handles.uibuttongroup1.SelectedObject.String
    case 'B' , SITE.SGS = 1;
    case 'C' , SITE.SGS = 2;
    case 'D' , SITE.SGS = 3;
end

varargout{1} = SUB;
varargout{2} = SC;
varargout{3} = SITE;


close(gcf)

function figure1_CloseRequestFcn(hObject, eventdata, handles) %#ok<*INUSD>
if isequal(get(hObject,'waitstatus'),'waiting')
    uiresume(hObject);
else
    delete(hObject);
end
% delete(hObject);

% ------------ subduction EQ ----------------------------------------------
function popupmenu1_Callback(hObject, eventdata, handles)

function edit1_Callback(hObject, eventdata, handles) %#ok<*DEFNU>

function edit2_Callback(hObject, eventdata, handles)

function edit3_Callback(hObject, eventdata, handles)

function edit4_Callback(hObject, eventdata, handles)

function edit5_Callback(hObject, eventdata, handles)

function edit6_Callback(hObject, eventdata, handles)


% ------------ subduction EQ ----------------------------------------------

function popupmenu2_Callback(hObject, eventdata, handles)

function edit11_Callback(hObject, eventdata, handles)

function edit12_Callback(hObject, eventdata, handles)

function edit13_Callback(hObject, eventdata, handles)

function edit14_Callback(hObject, eventdata, handles)

function edit15_Callback(hObject, eventdata, handles)

function edit16_Callback(hObject, eventdata, handles)

function edit17_Callback(hObject, eventdata, handles)

function edit18_Callback(hObject, eventdata, handles)

function edit19_Callback(hObject, eventdata, handles)

function edit20_Callback(hObject, eventdata, handles)

function edit21_Callback(hObject, eventdata, handles)

% ------------- site class ------------------------------------
function VS30_Callback(hObject, eventdata, handles)

function f0_Callback(hObject, eventdata, handles)

function Z10_Callback(hObject, eventdata, handles)

function Z25_Callback(hObject, eventdata, handles)

function uibuttongroup3_SelectionChangedFcn(hObject, eventdata, handles)

function popupmenu3_Callback(hObject, eventdata, handles)
