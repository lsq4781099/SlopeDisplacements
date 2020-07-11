function varargout = manageLayers(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @manageLayers_OpeningFcn, ...
    'gui_OutputFcn',  @manageLayers_OutputFcn, ...
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

function manageLayers_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>

handles.up.CData      = double(imread('up.jpg'))/255;
handles.down.CData    = double(imread('down.jpg'))/255;
handles.delete.CData  = double(imread('selection_delete.jpg'))/255;

handles.layer = varargin{1};
fld           = fields(handles.layer);
handles.layerpop.String =fld;
handles.listbox1.String = handles.layer.(fld{1});

guidata(hObject, handles);
uiwait(handles.figure1);

function varargout = manageLayers_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.layer;
delete(handles.figure1)

function listbox1_Callback(hObject, eventdata, handles) %#ok<*DEFNU>

function listbox1_ButtonDownFcn(hObject, eventdata, handles)

[filename,pathname,FILTERINDEX]=uigetfile({'*.mat','Source (*.mat)'},'Select a File','MultiSelect','on');

if FILTERINDEX==0
    return
end
if ischar(filename)
    filename = {filename};
end

if strcmp(handles.listbox1.String,' ')
    handles.listbox1.String = cell(1,length(filename));
    for i=1:length(filename)
        handles.listbox1.String{i}=fullfile(pathname,filename{i});
    end
else
    newString = cell(length(filename),1);
    for i=1:length(filename)
        newString{i}=fullfile(pathname,filename{i});
    end
    handles.listbox1.String=[handles.listbox1.String(1:end-1);newString;handles.listbox1.String(end)];
end

fld = handles.layerpop.String{handles.layerpop.Value};
handles.layer.(fld) = handles.listbox1.String;
guidata(hObject,handles)

function figure1_CloseRequestFcn(hObject, eventdata, handles)
if isequal(get(hObject,'waitstatus'),'waiting')
    uiresume(hObject);
else
    delete(hObject);
end
% delete(hObject);

function up_Callback(hObject, eventdata, handles)

Nflds = length(handles.listbox1.String);
if Nflds>=3 && handles.listbox1.Value<Nflds
    move = handles.listbox1.Value;
    str  = handles.listbox1.String;
    list = 1:length(str);
    if move>1
        list(move)=move-1;
        list(move-1)=move;
        handles.listbox1.String=str(list);
        handles.listbox1.Value=handles.listbox1.Value-1;
    end
    guidata(hObject,handles)
end

function down_Callback(hObject, eventdata, handles)
Nflds = length(handles.listbox1.String);
if Nflds>=3 && handles.listbox1.Value<Nflds-1
    move = handles.listbox1.Value;
    str  = handles.listbox1.String;
    list = 1:length(str);
    if move<length(str)
        list(move)=move+1;
        list(move+1)=move;
        handles.listbox1.String=str(list);
        handles.listbox1.Value=handles.listbox1.Value+1;
    end
    guidata(hObject,handles)
end

function delete_Callback(hObject, eventdata, handles)
Nflds = length(handles.listbox1.String);
if Nflds>=2 && handles.listbox1.Value~=Nflds
    del  = handles.listbox1.Value;
    str  = handles.listbox1.String;
    list = 1:length(str);
    list(del)=[];
    handles.listbox1.Value=1;
    handles.listbox1.String=str(list);
    
    fld = handles.layerpop.String{handles.layerpop.Value};
    handles.layer.(fld) = str(list);
    
    if isempty(str(list))
        handles.listbox1.String=' ';
        handles.listbox1.Value=1;
    end
    guidata(hObject,handles)
end

function layerpop_Callback(hObject, eventdata, handles)

fld = fields(handles.layer);
val = handles.layerpop.Value;
handles.listbox1.Value  = 1;
handles.listbox1.String = handles.layer.(fld{val});
