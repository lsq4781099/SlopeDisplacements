function varargout = VPSHA(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @VPSHA_OpeningFcn, ...
    'gui_OutputFcn',  @VPSHA_OutputFcn, ...
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

function VPSHA_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>
handles.ExitButton.CData = double(imread('exit.jpg'))/255;
handles.gridbutton.CData = double(imread('Grid.jpg'))/255;
handles.undock.CData      = double(imread('Undock.jpg'))/255;
handles.sys    = varargin{1};
handles.opt    = varargin{2};
handles.h      = varargin{3};

% ------------ REMOVES PCE MODELS (AT LEAST FOR NOW) ----------------------
isREG = handles.sys.isREG;
[~,B]=setdiff(handles.sys.branch(:,2),isREG);
if ~isempty(B)
    handles.sys.branch(B,:)=[];
    handles.sys.weight(B,:)=[];
    handles.sys.weight(:,5)=handles.sys.weight(:,5)/sum(handles.sys.weight(:,5));
    warndlg('PCE Models removed from logic tree. Logic Tree weights were normalized')
    uiwait
end

% -------------------------------------------------------------------------
Nmodels = size(handles.sys.branch,1);
handles.popbranch.String = compose('Branch %i',1:Nmodels);
handles.popsource.String = [{'all sources'};handles.sys.labelG{1}];
handles.popsite.String   = handles.h.id;

if numel(handles.opt.IM)==1
    answer = inputdlg('Specify 2nd intensity measure','VPSHA',1,{'SA(T=1)'});
    if isempty(answer)
        close(handles.figure1)
        return
    end
    handles.opt.IM(2)=str2IM(answer);
end

handles.popIM1.String    = IM2str(handles.opt.IM); handles.popIM1.Value=1;
handles.popIM2.String    = IM2str(handles.opt.IM); handles.popIM2.Value=2;

xlabel(handles.ax1,addIMunits(handles.opt.IM(1)),'fontsize',10)
ylabel(handles.ax1,addIMunits(handles.opt.IM(end)),'fontsize',10)
caxis([0 1])
C=colorbar('peer',handles.ax1,'location','eastoutside','tag','cbar');
set(get(C,'Title'),'String','MRE')

methods = pshatoolbox_methods(4);
handles.popcorrelation.String = {'Manual',methods(2:end).label};


akZoomTitle(handles.ax1)
guidata(hObject, handles);
% uiwait(handles.figure1);

function varargout = VPSHA_OutputFcn(hObject, eventdata, handles)
varargout{1} = [];

function popIM1_Callback(hObject, eventdata, handles) %#ok<*DEFNU>
methods = pshatoolbox_methods(4);
[~,ind] = intersect({methods.label},handles.popcorrelation.String{handles.popcorrelation.Value});

if isempty(ind)
    return
end

fun = methods(ind).func;
param.opp = 0; %
param.mechanism = 'interface';
param.direction = 'horizontal';
param.residual  = 'phi';
param.M         = 6;

IM1 = handles.opt.IM(handles.popIM1.Value);
IM2 = handles.opt.IM(handles.popIM2.Value);
rho = fun(IM1,IM2,param);
handles.rho.String = num2str(rho);

guidata(hObject,handles);

function popIM2_Callback(hObject, eventdata, handles)
methods = pshatoolbox_methods(4);
[~,ind] = intersect({methods.label},handles.popcorrelation.String{handles.popcorrelation.Value});

if isempty(ind)
    return
end

fun = methods(ind).func;
param.opp = 0; %
param.mechanism = 'interface';
param.direction = 'horizontal';
param.residual  = 'phi';
param.M         = 6;

IM1 = handles.opt.IM(handles.popIM1.Value);
IM2 = handles.opt.IM(handles.popIM2.Value);
rho = fun(IM1,IM2,param);
handles.rho.String = num2str(rho);

guidata(hObject,handles);

function run2Bazurro_Callback(hObject, eventdata, handles)

model_ptr  = handles.popbranch.Value;
source_ptr = handles.popsource.Value;
site_ptr   = handles.popsite.Value;
IM_ptr     = [handles.popIM1.Value,handles.popIM2.Value];
branch     = handles.sys.branch;

if source_ptr ==1
    Nsource    = sum(handles.sys.Nsrc(:,branch(model_ptr,1)));
    source_ptr = 1:Nsource;
else
    source_ptr = source_ptr-1;
end

ncol    = size(handles.opt.im,2);
if ncol==1
    im      = repmat(handles.opt.im,1,2);
else
    im      = handles.opt.im(:,IM_ptr);
end

IM      = handles.opt.IM(IM_ptr);
h.id    = handles.h.id(site_ptr,:);
h.p     = handles.h.p(site_ptr,:);
h.param = handles.h.param;
h.value = handles.h.value(site_ptr,:);

opt     = handles.opt;
sources = buildmodelin(handles.sys,branch(model_ptr,:),handles.opt);
sources = sources(source_ptr);
Nsource = length(source_ptr);
rho     = str2double(handles.rho.String);
[handles.lambda,handles.MRE,handles.MRD] = runhazardV1(im,IM,h,opt,sources,Nsource,rho);

handles.im = im;
handles.IM = IM;
plotVPSHA(handles)
guidata(hObject,handles)

function popsite_Callback(hObject, eventdata, handles) %#ok<*INUSD>

function popbranch_Callback(hObject, eventdata, handles)

function popsource_Callback(hObject, eventdata, handles)

function rho_Callback(hObject, eventdata, handles)

function popcorrelation_Callback(hObject, eventdata, handles)

methods = pshatoolbox_methods(4);
[~,ind] = intersect({methods.label},handles.popcorrelation.String{handles.popcorrelation.Value});

if isempty(ind)
    return
end

fun = methods(ind).func;
param.opp = 0; %
param.mechanism = 'interface';
param.direction = 'horizontal';
param.residual  = 'phi';
param.M         = 6;

IM1 = handles.opt.IM(handles.popIM1.Value);
IM2 = handles.opt.IM(handles.popIM2.Value);
rho = fun(IM1,IM2,param);
handles.rho.String = num2str(rho);

guidata(hObject,handles);

function uibuttongroup1_SelectionChangedFcn(hObject, eventdata, handles)
plotVPSHA(handles)

function figure1_CloseRequestFcn(hObject, eventdata, handles)
delete(hObject);

function File_Callback(hObject, eventdata, handles)

function Exit_Callback(hObject, eventdata, handles)
close(handles.figure1)

function Advanced_Callback(hObject, eventdata, handles)

function ptype1_Callback(hObject, eventdata, handles)
switch hObject.Checked
    case 'off'
        hObject.Checked='on';
        handles.ptype2.Checked='off';
        handles.ptype3.Checked='off';
        plotVPSHA(handles)
end

function ptype2_Callback(hObject, eventdata, handles)
switch hObject.Checked
    case 'off'
        hObject.Checked='on';
        handles.ptype1.Checked='off';
        handles.ptype3.Checked='off';
        plotVPSHA(handles)
end

function ptype3_Callback(hObject, eventdata, handles)
switch hObject.Checked
    case 'off'
        hObject.Checked='on';
        handles.ptype1.Checked='off';
        handles.ptype2.Checked='off';
        plotVPSHA(handles)
end

function ExitButton_Callback(hObject, eventdata, handles)
close(handles.figure1);

function gridbutton_Callback(hObject, eventdata, handles)

ax = handles.ax1;
switch [ax.XGrid,ax.XMinorGrid]
    case 'onon'  , ax.XGrid='on';  ax.YGrid='on';  ax.ZGrid='on';  ax.XMinorGrid='off'; ax.YMinorGrid='off'; ax.ZMinorGrid='off';
    case 'onoff' , ax.XGrid='off'; ax.YGrid='off'; ax.ZGrid='off'; ax.XMinorGrid='on';  ax.YMinorGrid='on';  ax.ZMinorGrid='on';
    case 'offon' , ax.XGrid='off'; ax.YGrid='off'; ax.ZGrid='off'; ax.XMinorGrid='off'; ax.YMinorGrid='off'; ax.ZMinorGrid='off';
    case 'offoff', ax.XGrid='on';  ax.YGrid='on';  ax.ZGrid='on';  ax.XMinorGrid='on';  ax.YMinorGrid='on';  ax.ZMinorGrid='on';
end

function undock_Callback(hObject, eventdata, handles)
figure2clipboard_uimenu(hObject, eventdata,handles.ax1)

function scalarhazard_Callback(hObject, eventdata, handles)

ch = findall(handles.ax1,'tag','pat2');
switch hObject.Checked
    case 'on' , hObject.Checked='off'; if ~isempty(ch),set(ch,'Visible','off'); end
    case 'off', hObject.Checked='on';  if ~isempty(ch),set(ch,'Visible','on');  end
end
guidata(hObject,handles)
