function varargout = SeismicHazard(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @SeismicHazard_OpeningFcn, ...
    'gui_OutputFcn',  @SeismicHazard_OutputFcn, ...
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

function SeismicHazard_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>
home
fprintf(' --------- Seismic Hazard Toolbox ---------\n')
load pshabuttons c1 c2 c3 c4 c5 c6 c7 c8 c9 c10 c11 c12 c13
handles.form1                      = c1;
handles.form2                      = c2;
handles.engine.CData               = c3;
handles.switchmode.CData           = c4;
handles.addLeg.CData               = c5;
handles.ax2Limits.CData            = c6;
handles.Exit_button.CData          = c7;
handles.Distance_button.CData      = c8;
handles.po_refresh_GE.CData        = c9;
handles.RefreshButton.CData        = c10;
handles.OpenRef.CData              = c11;
handles.ColorSecondaryLines.CData  = c12;
handles.restore_axis.CData         = c13;
handles.GoogleEarthOpt             = GEOptions('default');
initializePSHA

if nargin==4
    fname = varargin{1};
    [handles.sys,handles.opt,handles.h]=psha_updatemodel([],fname);
    house_keeping1
    handles.launch_PSDA.Enable='off';
    plot_sites_PSHA
end
akZoom(handles.ax1)
guidata(hObject, handles);
% uiwait(handles.fig);

function varargout = SeismicHazard_OutputFcn(~, ~, ~)
varargout{1}=[];
% delete(handles.fig)

function fig_CloseRequestFcn(hObject, ~, ~)
% if isequal(get(hObject,'waitstatus'),'waiting')
%     uiresume(hObject);
% else
%     delete(hObject);
% end
delete(hObject);

% -----FILE MENU ----------------------------------------------------------
function File_Callback(~, ~, ~) %#ok<*DEFNU>

function Reset_Callback(hObject, ~, handles)
initializePSHA;
guidata(hObject, handles);

function Clear_Analysis_Callback(hObject, eventdata, handles)
delete_analysis_PSHA
plot_sites_PSHA
guidata(hObject, handles);

function InspectInputFile_Callback(hObject, eventdata, handles)
if ~isempty(handles.sys)
    if ispc,winopen(handles.sys.filename);end
end

function LoadTXT_Callback(~, ~, ~)

function DefaultSeismicity_Callback(hObject, ~, handles)
anw = listdlg('PromptString','Default Seismicity Models:',...
    'SelectionMode','single','ListSize',[200 300],...
    'ListString',{...
    'Chile (Martin 1990)'...
    'Chile (Poulos et al., 2018)'...
    'Peru (SENCICO, 2016)',...
    'Mexico (2019)',...
    'Ciudad Universitaria, Mexico (2019)',...
    'Ecuador',...
    'Colombia',...
    'USGS (NSHM 2008 Dynamic)',...
    'USGS (NSHM 2014 Dynamic)'});
if isempty(anw)
    return
end
initializePSHA;
switch anw
    case 1,  [handles.sys,handles.opt,handles.h]=psha_updatemodel('Chile_Default'); house_keeping1;
    case 2,  [handles.sys,handles.opt,handles.h]=psha_updatemodel('Chile_Poulos');  house_keeping1;
    case 3,  [handles.sys,handles.opt,handles.h]=psha_updatemodel('Peru_Default');  house_keeping1;
    case 4,  [handles.sys,handles.opt,handles.h]=psha_updatemodel('Mexico_Default_not_CU'); house_keeping1;
    case 5,  [handles.sys,handles.opt,handles.h]=psha_updatemodel('MexicoCity_CU_PreviousTectonic_2018_f2'); house_keeping1;
    case 6,  [handles.sys,handles.opt,handles.h]=psha_updatemodel('Ecuador_Default'); house_keeping1;
    case 7,  [handles.sys,handles.opt,handles.h]=psha_updatemodel('Colombia_Default'); house_keeping1;
    case 8,  handles=usgs_updatemodel(handles,'USGS_NHSM_2008');
    case 9,  handles=usgs_updatemodel(handles,'USGS_NHSM_2014');
end
plot_sites_PSHA
guidata(hObject,handles)

function TestModels_Callback(hObject, eventdata, handles)
W = what('platform_testmodels');
D = dir(W.path);
D = D(contains({D.name},'.txt'));
anw = listdlg('PromptString','Test Models:',...
    'SelectionMode','single','ListString',{D.name});
if isempty(anw)
    return
end
initializePSHA;
[handles.sys,handles.opt,handles.h]=psha_updatemodel([],fullfile(W.path,D(anw).name));
house_keeping1
plot_sites_PSHA

guidata(hObject, handles);

function ValidationPEER2018_Callback(hObject,envetdata,handles)
anw = listdlg('PromptString','PEER 2018 Validation Tests:',...
    'SelectionMode','single','ListString',{'Test 1.1','Test 1.2','Test 1.3','Test 1.4','Test 1.5',...
    'Test 1.6','Test 1.7','Test 1.8a','Test 1.8b','Test 1.8c','Test 1.10','Test 1.11',...
    'Test 2.1','Test 2.2a','Test 2.2b','Test 2.2c','Test 2.2d',...
    'Test 2.3a','Test 2.3b','Test 2.3c','Test 2.3d'});
if isempty(anw)
    return
end
initializePSHA;
switch anw
    case 1,  [handles.sys,handles.opt,handles.h]=psha_updatemodel('PEER2018_Set1_Test1_1');
    case 2,  [handles.sys,handles.opt,handles.h]=psha_updatemodel('PEER2018_Set1_Test1_2');
    case 3,  [handles.sys,handles.opt,handles.h]=psha_updatemodel('PEER2018_Set1_Test1_3');
    case 4,  [handles.sys,handles.opt,handles.h]=psha_updatemodel('PEER2018_Set1_Test1_4');
    case 5,  [handles.sys,handles.opt,handles.h]=psha_updatemodel('PEER2018_Set1_Test1_5');
    case 6,  [handles.sys,handles.opt,handles.h]=psha_updatemodel('PEER2018_Set1_Test1_6');
    case 7,  [handles.sys,handles.opt,handles.h]=psha_updatemodel('PEER2018_Set1_Test1_7');
    case 8,  [handles.sys,handles.opt,handles.h]=psha_updatemodel('PEER2018_Set1_Test1_8a');
    case 9,  [handles.sys,handles.opt,handles.h]=psha_updatemodel('PEER2018_Set1_Test1_8b');
    case 10, [handles.sys,handles.opt,handles.h]=psha_updatemodel('PEER2018_Set1_Test1_8c');
    case 11, [handles.sys,handles.opt,handles.h]=psha_updatemodel('PEER2018_Set1_Test1_10');
    case 12, [handles.sys,handles.opt,handles.h]=psha_updatemodel('PEER2018_Set1_Test1_11');
    case 13, [handles.sys,handles.opt,handles.h]=psha_updatemodel('PEER2018_Set2_Test2_1');
    case 14, [handles.sys,handles.opt,handles.h]=psha_updatemodel('PEER2018_Set2_Test2_2a');
    case 15, [handles.sys,handles.opt,handles.h]=psha_updatemodel('PEER2018_Set2_Test2_2b');
    case 16, [handles.sys,handles.opt,handles.h]=psha_updatemodel('PEER2018_Set2_Test2_2c');
    case 17, [handles.sys,handles.opt,handles.h]=psha_updatemodel('PEER2018_Set2_Test2_2d');
    case 18, [handles.sys,handles.opt,handles.h]=psha_updatemodel('PEER2018_Set2_Test2_3a');
    case 19, [handles.sys,handles.opt,handles.h]=psha_updatemodel('PEER2018_Set2_Test2_3b');
    case 20, [handles.sys,handles.opt,handles.h]=psha_updatemodel('PEER2018_Set2_Test2_3c');
    case 21, [handles.sys,handles.opt,handles.h]=psha_updatemodel('PEER2018_Set2_Test2_3d');
end
house_keeping1
plot_sites_PSHA
handles.pdffile = handles.sys.filename;
guidata(hObject,handles)

function LoadCustom_Callback(hObject, eventdata, handles)

if isfield(handles,'defaultpath_others')
    [filename,pathname,FILTERINDEX]=uigetfile(fullfile(handles.defaultpath_others,'*.txt'),'select file');
else
    [filename,pathname,FILTERINDEX]=uigetfile(fullfile(pwd,'*.txt'),'select file');
end

if FILTERINDEX==0,return;end
initializePSHA;
[handles.sys,handles.opt,handles.h]=psha_updatemodel([],[pathname,filename]);
house_keeping1
plot_sites_PSHA
handles.defaultpath_others=pathname;
guidata(hObject, handles);

function Export_Results_Callback(~,~,~)

function ExportHazard_Callback(hObject, eventdata, handles)
[FileName,PathName,FILTERINDEX] =  uiputfile('*.out','Save Hazard Analysis As');
if FILTERINDEX==0
    return
end
WEIGHT = handles.sys.weight(:,5);
haz2txt(FileName,PathName,WEIGHT,handles.MRE,handles.opt,handles.h)

function ExportIS_Callback(~, ~, handles) %#ok<*INUSD>
[FileName,PathName,FILTERINDEX] =  uiputfile('*.bin','Save Hazard Analysis As'); %#ok<*ASGLU>
FILE = fullfile(PathName,FileName);
if FILTERINDEX==0
    return
end
etype='IS';
dsha_bin

function ExportKM_Callback(hObject, eventdata, handles)
[FileName,PathName,FILTERINDEX] =  uiputfile('*.bin','Save Hazard Analysis As');
FILE = fullfile(PathName,FileName);
if FILTERINDEX==0
    return
end
etype='KM';
dsha_bin

function ExportGE_Image_Callback(hObject, eventdata, handles)

D=what('DEFAULT_MATFILES');
if isempty(D)
    [filename,pathname,FILTERINDEX] =  uiputfile('*.mat','Save Google Earth Image');
else
    [filename,pathname,FILTERINDEX] =  uiputfile([D.path,'\*.mat'],'Save Google Earth Image');
end
if FILTERINDEX==0
    return
end
ch=findall(handles.ax1,'tag','gmap');
if ~isempty(ch)
    xx = ch.XData; yy = ch.YData; cc = ch.CData; %#ok<*NASGU>
    XLIM=handles.ax1.XLim;
    YLIM=handles.ax1.YLim;
    save([pathname,filename],'xx','yy','cc','XLIM','YLIM');
else
    warndlg('No Background Image Found')
end

function Undockax_Callback(~,~,~)

function ExportHazardCurves_Callback(hObject, eventdata, handles)
figure2clipboard_uimenu(hObject, eventdata,handles.ax2)

function ExportMap_Callback(hObject, eventdata, handles)
figure2clipboard_uimenu(hObject, eventdata,handles.ax1)

function Exit_Callback(~,~,~)
close(gcf)

function Exit_button_Callback(hObject, eventdata,~)
close(gcf)

% ------EDIT MENU --------------------------------------------------------
function Edit_Callback(~, ~, ~)

function Update_Global_Parameters_Callback(hObject, ~, handles)

old_opt     = handles.opt;
switch handles.sys.filename
    case {'USGS_NHSM_2008','USGS_NHSM_2014'}
        new_opt     = GlobalParamUSGS(handles.opt);
        handles.opt = new_opt;
        handles.IM_select.String=handles.opt.IM;
        if ~structcmp(new_opt,old_opt)
            delete_analysis_PSHA;
        end
    otherwise
        new_opt = GlobalParam(handles.opt);
        
        if strcmp(new_opt.SourceDeagg,'off')
            handles.HazOptions.sbh(2:3)=0;
        end
        A=unique(handles.sys.gmmptr(:));
        Tmax = [];
        for i=1:length(A)
            T    = handles.sys.gmmlib(A(i)).T;
            Tmax = [Tmax;max(T)]; %#ok<AGROW>
        end
        Tmax      = min(Tmax);
        Tdeleted  = new_opt.IM(new_opt.IM>Tmax);
        new_opt.IM(new_opt.IM>Tmax) = [];
        
        if ~isempty(Tdeleted)
            h=warndlg({['Analysis for T = ',mat2str(Tdeleted'),'  not available'],...
                'check limits of GMMs'});
            uiwait(h);
        end
        
        handles.opt  = new_opt;
        handles.IM_select.Value  = 1;
        handles.IM_select.String = IM2str(handles.opt.IM);
        handles.pop_field.Value  = 1;
        handles.pop_field.String = IM2str([handles.opt.IM1;handles.opt.IM2(:)]);
        
        if ~structcmp(new_opt,old_opt)
            delete_analysis_PSHA;
        end
        if strcmp(handles.opt.Clusters{1},'on')
            handles.po_clusters.Enable='on';
            [handles.idx,handles.hc] = compute_clusters(handles.opt,handles.h);
        else
            handles.idx=[];
            handles.hc=createObj('site');
            handles.po_clusters.Enable='off';
        end
        plot_sites_PSHA
        default_maps(handles.fig,handles.opt);
end
guidata(hObject, handles);

function LogicTree_Callback(hObject, eventdata, handles)
if ~~isempty(handles.sys), return;end
handles.sys.branch=PSHA_Logic_Tree(handles.sys.branch,handles.sys.weight,handles.sys.gmmid);
guidata(hObject,handles)

function SourceGeometry_Callback(hObject, eventdata, handles)
SourceGeometry(handles);
guidata(hObject, handles);

function MagnitudeRecurrence_Callback(hObject, eventdata, handles)
MR_explorer(...
    handles.sys.mrr1,...
    handles.sys.mrr2,...
    handles.sys.mrr3,...
    handles.sys.mrr4,...
    handles.sys.mrr5,...
    handles.sys.Nmrr,...
    handles.opt)

function Seismicity_Callback(hObject, eventdata, handles)
GMM_explorer(handles)

function Sites_Callback(hObject, eventdata, handles)

handin.h            = handles.h;
handin.ax           = handles.ax1;
handin.opt          = handles.opt;
handin.layer        = handles.sys.layer;
[h_new,layer,index] = SelectLocations(handin);
equalh              = structcmp(handin.h,h_new);

switch handles.opt.Clusters{1}
    case 'on' , [handles.idx,handles.hc] = compute_clusters(handles.opt,handles.h);
    case 'off'
        handles.idx=[];
        handles.hc =createObj('site');
end

if equalh
    return
end

delete_analysis_PSHA;
handles.site_selection = 1:length(h_new.id);
Nsites = size(h_new.p,1);
handles.h = h_new;
handles.sys.layer = layer;
handles.index = index;
handles.site_menu.String = handles.h.id;
handles.site_menu.Value  = 1;
handles.site_menu_psda.String = handles.h.id;
handles.site_menu_psda.Value  = 1;
plot_sites_PSHA

guidata(hObject, handles);

function GEmaps_Callback(hObject, eventdata, handles)

if ~exist('api_Key.mat','file')
    warndlg('You must use an API key to authenticate each request to Google Maps Platform APIs. For additional information, please refer to http://g.co/dev/maps-no-account')
    return
end

handles.GoogleEarthOpt=GEOptions(handles.GoogleEarthOpt);
guidata(hObject, handles);

% -----TOOLS MENU --------------------------------------------------------
function SeismicHazardTools_Callback(~,~,~)

function launch_UHS_Callback(hObject, eventdata, handles)
if ~isempty(handles.sys) && ~isempty(handles.sys.isREG)
    UHS(handles.sys,handles.opt,handles.h)
end

function launch_CMS_Callback(hObject, eventdata, handles)
if ~isempty(handles.sys) && ~isempty(handles.sys.isREG)
    CMS(handles.sys,handles.opt,handles.h)
end

function launch_CSS_Callback(hObject, eventdata, handles)
if ~isempty(handles.sys) && ~isempty(handles.sys.isREG)
    CSS(handles.sys,handles.opt,handles.h)
end

function launch_CSSIa_Callback(hObject, eventdata, handles)
if ~isempty(handles.sys) && ~isempty(handles.sys.isREG)
    %     CSSAI(handles.sys,handles.opt,handles.h)
end

function launch_DEAGG_MR_Callback(hObject, eventdata, handles)
switch handles.sys.filename
    case {'USGS_NHSM_2008','USGS_NHSM_2014'}
        % DeaggregationUSGS(handles); % UNDER MAINTENANCE
    otherwise
        if ~isempty(handles.sys) && ~isempty(handles.sys.isREG)
            DeaggregationMR(handles.sys,handles.opt,handles.h);
        end
end
guidata(hObject,handles)

function launch_DEAGG_MRe_Callback(hObject, eventdata, handles)
switch handles.sys.filename
    case {'USGS_NHSM_2008','USGS_NHSM_2014'} % not available
    otherwise
        if ~isempty(handles.sys) && ~isempty(handles.sys.isREG)
            DeaggregationMRe(handles.sys,handles.opt,handles.h);
        end
end
guidata(hObject,handles)

function launch_PSDA_Callback(hObject, eventdata, handles)
PSDA2(handles.sys,handles.opt,handles.h)

function launch_PVSHA_Callback(hObject, eventdata, handles)
if ~isempty(handles.sys) && ~isempty(handles.sys.isREG)
    VPSHA(handles.sys,handles.opt,handles.h)
end

function launch_PCE_Callback(hObject, eventdata, handles)
if ~isempty(handles.sys.isPCE)
    CGMM(handles.sys,handles.opt,handles.h,handles.HazOptions)
end

function launch_GRD_Callback(hObject, eventdata, handles)
if ~isempty(handles.sys) && ~isempty(handles.sys.isREG)
    GRD(handles.sys,handles.opt,handles.h,handles.ax1)
end

% -----PSHA MENU --------------------------------------------------------
function PSHAsetup_Callback(hObject, eventdata, handles)
handles.site_selection=PSHARunOptions(handles);
guidata(hObject,handles)

function PSHAMenu_Callback(~,~,~)

function site_menu_Callback(hObject, eventdata, handles)
plot_hazard_PSHA(...
    handles.fig,...
    handles.HazOptions,...
    handles.opt,...
    handles.sys.isPCE,...
    handles.MRE,...
    handles.MREPCE,...
    handles.sys.validation,...
    handles.sys.weight(:,5),...
    handles.sys.labelG,...
    handles.sys.mech,...
    handles.idx,1);

guidata(hObject, handles);

function IM_select_Callback(hObject, eventdata, handles)
plot_hazard_PSHA(...
    handles.fig,...
    handles.HazOptions,...
    handles.opt,...
    handles.sys.isPCE,...
    handles.MRE,...
    handles.MREPCE,...
    handles.sys.validation,...
    handles.sys.weight(:,5),...
    handles.sys.labelG,...
    handles.sys.mech,...
    handles.idx);

IM_ptr = handles.IM_select.Value;
if ismember(handles.HazOptions.sbh(1),handles.sys.isREG)
    handles.site_colors = compute_v(...
        handles.fig,...
        handles.opt,...
        handles.HazOptions,...
        nansum(handles.MRE(:,:,IM_ptr,:,:),4),...
        handles.sys.weight(:,5));
    plot_hazardmap_PSHA(...
        handles.fig,...
        handles.site_colors,...
        handles.HazOptions,...
        handles.h,...
        handles.idx);
end
guidata(hObject, handles);

function runMRE_Callback(hObject, eventdata, handles)

opt    = handles.opt;
slist1 = handles.site_selection;
slist2 = 1:max(handles.idx);

switch handles.sys.filename
    case 'USGS_NHSM_2008', handles=runhazUSGS(handles);
    case 'USGS_NHSM_2014', handles=runhazUSGS(handles);
    otherwise
        switch handles.opt.Clusters{1}
            case 'off', [handles.MRE,handles.MREPCE] = runlogictree1(handles.sys,opt,handles.h ,slist1);
            case 'on' , [handles.MRE,handles.MREPCE] = runlogictree1(handles.sys,opt,handles.hc,slist2);
        end
        
        handles.site_menu.Value      = handles.site_selection(1);
        plot_hazard_PSHA(...
            handles.fig,...
            handles.HazOptions,...
            handles.opt,...
            handles.sys.isPCE,...
            handles.MRE,...
            handles.MREPCE,...
            handles.sys.validation,...
            handles.sys.weight(:,5),...
            handles.sys.labelG,...
            handles.sys.mech,...
            handles.idx);
        
        IM_ptr = handles.IM_select.Value;
        handles.site_colors = compute_v(...
            handles.fig,...
            handles.opt,...
            handles.HazOptions,...
            nansum(handles.MRE(:,:,IM_ptr,:,:),4),...
            handles.sys.weight(:,5));
        plot_hazardmap_PSHA(...
            handles.fig,...
            handles.site_colors,...
            handles.HazOptions,...
            handles.h,...
            handles.idx);
end
handles.platformMode = 1;
handles.PSHApannel.Visible='on';
handles.DSHApannel.Visible='off';
handles.switchmode.Enable=enableswitchmode(~isempty(handles.MRE),~isempty(handles.scenarios));
handles.switchmode.CData=handles.form1;
handles.switchmode.TooltipString='switch to PSDA mode';
handles.ax2.XScale='log';
handles.ax2.YScale='log';

if ispc && isfield(handles.sys,'validation') && ~isempty(handles.sys.validation)
    comparePEER(handles.sys.validation,handles.MRE,opt.im,handles.sys.filename)
end

if ispc
    pdffile = strrep(handles.sys.filename,'.txt','.pdf');
    if exist(pdffile,'file')
        handles.pdffile = pdffile;
        handles.OpenRef.Visible='on';
    else
        handles.OpenRef.Visible='off';
    end
end

guidata(hObject, handles);

% ---------------- DSHAMENU MENU----------------------------
function DSHAMenu_Callback(~,~,~)

function CreateScenarios_Callback(hObject, eventdata, handles)

[handles.scenarios,handles.source] = ScenarioBuilder2(handles.sys,handles.opt,handles.h,handles.scenarios,handles.source);
Ns = size(handles.scenarios,1);
handles.pop_scenario.Value   = 1;
handles.pop_scenario.String  = compose('S%i',1:Ns);
handles.pop_field.String     = IM2str([handles.opt.IM1(:);handles.opt.IM2(:)]);
guidata(hObject,handles)

function Run_DHSA_IS_Callback(hObject, eventdata, handles)

if isempty(handles.source) || isempty(handles.scenarios)
    handles.source    = buildDSHAmodel2(handles.sys,handles.opt,handles.h);
    handles.scenarios = sSCEN_list(handles.source);
    Ns = size(handles.scenarios,1);
    handles.pop_scenario.Value   = 1;
    handles.pop_scenario.String  = compose('S%i',1:Ns);
    handles.pop_field.String     = IM2str([handles.opt.IM1(:);handles.opt.IM2(:)]);
end

handles.pop_scenario.Value = 1;
[handles.mulogIM,handles.L] = run_DSHA(handles.source,handles.scenarios,handles.h,handles.opt);
plot_DSHA

handles.platformMode = 2;
handles.PSHApannel.Visible='off';
handles.DSHApannel.Visible='on';
handles.switchmode.Enable=enableswitchmode(~isempty(handles.MRE),~isempty(handles.scenarios));
handles.switchmode.CData=handles.form2;
handles.switchmode.TooltipString='switch to PSHA mode';
dsha_lambda2(handles);

guidata(hObject,handles)

function Run_DSHA_KM_Callback(hObject, eventdata, handles)
answer = inputdlg({'Number of Clusters:','Replicate'},'k-means',[1,35],{'200','1'});
if isempty(answer)
    return
end
handles.optkm = str2double(answer);
[handles.krate,handles.kY] = dsha_kmeans(handles,handles.optkm);
dsha_lambda2(handles);
guidata(hObject,handles)

function pop_scenario_Callback(hObject, eventdata, handles)
plot_DSHA

function pop_field_Callback(hObject, eventdata, handles)
plot_DSHA
dsha_lambda2(handles);
guidata(hObject,handles)

function RefreshButton_Callback(hObject, eventdata, handles)
plot_DSHA
dsha_lambda2(handles);
guidata(hObject,handles)

function show_RA_Callback(hObject, eventdata, handles)
ch = findall(handles.ax1,'tag','scenario');
switch hObject.Value
    case 1,set(ch,'Visible','on');
    case 2,set(ch,'Visible','off');
end
guidata(hObject,handles)

function PSDA_display_mode_Callback(hObject, eventdata, handles)
plot_DSHA
guidata(hObject,handles)

function NumSim_Callback(hObject, eventdata, handles)
dsha_lambda2(handles);
guidata(hObject,handles)

function shakecorrmode_Callback(hObject, eventdata, handles)
opt2 = handles.opt;
switch hObject.Value
    case 1 % default
    case 2, opt2.Spectral=@none_spectral;
    case 3, opt2.Spatial=@none_spatial;
    case 4, opt2.Spatial=@none_spatial; opt2.Spectral=@none_spectral;
end
[handles.mulogIM,handles.L] = run_DSHA(handles.source,handles.scenarios,handles.h,opt2);
plot_DSHA
dsha_lambda2(handles);
guidata(hObject,handles)

% -----DISPLAY OPTIONS ------------------------------------------------------
function po_grid_Callback(hObject, eventdata, handles)

val = hObject.Value;
switch val
    case 1, grid(handles.ax1,'on')
    case 0, grid(handles.ax1,'off')
end
guidata(hObject,handles);

function Boundary_check_Callback(hObject, eventdata, handles)
ch = findall(handles.ax1,'tag','shape1');
val = hObject.Value;
switch val
    case 1, set(ch,'visible','on');
    case 0, set(ch,'visible','off')
end
guidata(hObject,handles);

function po_sources_Callback(hObject, eventdata, handles)
i   = handles.po_region.Value;
tag_edge = sprintf('edge%g',handles.po_region.Value);
tag_mesh = sprintf('mesh%g',handles.po_region.Value);
switch hObject.Value
    case 1
        ch=findall(handles.ax1,'tag',tag_edge); set(ch,'visible','on');
    case 0
        ch=findall(handles.ax1,'tag',tag_edge); set(ch,'visible','off');
        ch=findall(handles.ax1,'tag',tag_mesh); set(ch,'visible','off');
        handles.po_sourcemesh.Value=0;
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

function po_sites_Callback(hObject, eventdata, handles)
ch=findall(handles.ax1,'tag','siteplot');
switch hObject.Value
    case 1, set(ch,'visible','on');
    case 0, set(ch,'visible','off');
end
guidata(hObject,handles);

function po_clusters_Callback(hObject, eventdata, handles)
ch=findall(handles.ax1,'tag','cluster');
switch hObject.Value
    case 1, set(ch,'visible','on');
    case 0, set(ch,'visible','off');
end

function po_contours_Callback(hObject, eventdata, handles)
if isempty(handles.MRE) && isempty(handles.scenarios)
    return
end
switch handles.po_contours.Value
    case 0
        for i=1:size(handles.h.t,1)
            ch = findall(handles.ax1,'tag',num2str(i));
            ch.Visible='off';
        end
        ch=findall(handles.ax1,'Tag','satext');
        if ~isempty(ch)
            ch.Visible='off';
        end
    case 1
        for i=1:size(handles.h.t,1)
            ch = findall(handles.ax1,'tag',num2str(i));
            ch.Visible='on';
        end
        ch=findall(handles.ax1,'Tag','satext');
        if ~isempty(ch)
            ch.Visible='on';
        end
end
v = handles.po_sites.Value+handles.po_contours.Value;
ch = findall(handles.fig,'type','ColorBar');
if ~isempty(ch)
    switch v
        case 0    ,ch.Visible='off';
        otherwise ,ch.Visible='on';
    end
end
guidata(hObject, handles);

function po_region_Callback(hObject, eventdata, handles)

ngeom = length(hObject.String);
for i=1:ngeom
    tag=sprintf('edge%g',i);ch=findall(handles.ax1,'tag',tag);set(ch,'Visible','off');
    tag=sprintf('mesh%g',i);ch=findall(handles.ax1,'tag',tag);set(ch,'Visible','off');
end


for j=1:ngeom
    tag=sprintf('sourcetag%g',j);
    ch = findall(handles.ax1,'tag',tag);
    set(ch,'Visible','off');
end

if handles.po_sources.Value
    tag=sprintf('edge%g',hObject.Value);
    ch=findall(handles.ax1,'tag',tag);
    set(ch,'Visible','on');
end

if handles.po_sourcemesh.Value
    tag=sprintf('mesh%g',hObject.Value);
    ch=findall(handles.ax1,'tag',tag);
    set(ch,'Visible','on');
end

if handles.SourceLabels.Value
    tag=sprintf('sourcetag%g',hObject.Value);
    ch = findall(handles.ax1,'tag',tag);
    set(ch,'Visible','on');
end

function po_googleearth_Callback(hObject, eventdata, handles)
ch=findall(handles.ax1,'tag','gmap');
if isempty(ch)
    return
end
switch hObject.Value
    case 0, ch.Visible='off';
    case 1, ch.Visible='on';
end

function po_refresh_GE_Callback(hObject, eventdata, handles)
if handles.opt.ellipsoid.Code==0
    return
end

if ~exist('api_Key.mat','file')
    warndlg('You must use an API key to authenticate each request to Google Maps Platform APIs. For additional information, please refer to http://g.co/dev/maps-no-account')
    return
end

a=1;
if ispc
    a=system('ping -n 1 www.google.com');
elseif isunix
    a=system('ping -c 1 www.google.com');
end
if a==1
    fprintf('Warning: No Internet Access found\n');
    return
end

opt=handles.GoogleEarthOpt;
caxis(handles.ax1);
gmap=plot_google_map(...
    'Axis',handles.ax1,...
    'Height',640,...
    'Width',640,...
    'Scale',2,...
    'MapType',opt.MapType,...
    'Alpha',opt.Alpha,...
    'ShowLabels',opt.ShowLabels,...
    'autoAxis',0,...
    'refresh',0,...
    'Color',opt.Color);

gmap.Tag='gmap';
ch = findall(handles.ax1,'Tag','gmap');
if length(ch)>1
    delete(ch(1))
    ch(1)=[];
end
try %#ok<TRYNC>
    uistack(gmap,'bottom');
    handles.ax1.Layer='top';
end
guidata(hObject,handles)

% ---------CUSTOM TOOLS----------------------------------------------------
function Distance_button_Callback(hObject, eventdata, handles)
ch1=findall(handles.ax1,'tag','patchselect');
ch2=findall(handles.ax1,'tag','patchtxt');
if isempty(ch1) && isempty(ch2)
    ch1=findall(handles.fig,'Style','pushbutton','Enable','on'); set(ch1,'Enable','inactive');
    ch2=findall(handles.fig,'type','uimenu','Enable','on'); set(ch2,'Enable','off');
    XYLIM1 = get(handles.ax1,{'xlim','ylim'});
    if handles.opt.ellipsoid.Code==0
        show_distanceECEF(handles.ax1,'line');
    else
        show_distance(handles.ax1,'line');
    end
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

function CreateAPI_Callback(~,~,~)

prompt = {'Paste API key:'};
title = 'API key';
dims = [1 50];
if exist('api_key.mat','file')
    load api_key.mat apiKey
    definput = {apiKey};
else
    definput = {''};
end
answer = inputdlg(prompt,title,dims,definput);
if isempty(answer)
    return
else
    apiKey = answer{1};
    D=what('platform_apiKey');
    filename = fullfile(D.path,'api_key.mat');
    save(filename,'apiKey');
end

function engine_Callback(hObject, eventdata, handles)

if isempty(handles.MRE)
    return
end
nmodels=size(handles.sys.branch,1);
mtype  = false(1,nmodels);
for i=1:length(handles.sys.isREG)
    mtype(i)=ismember(handles.sys.isREG(i),1:nmodels);
end
str = [compose('Branch %i',1:nmodels);num2cell(mtype)];
handles.HazOptions = HazardOptions(str,handles.HazOptions,handles.opt.SourceDeagg);

plot_hazard_PSHA(...
    handles.fig,...
    handles.HazOptions,...
    handles.opt,...
    handles.sys.isPCE,...
    handles.MRE,...
    handles.MREPCE,...
    handles.sys.validation,...
    handles.sys.weight(:,5),...
    handles.sys.labelG,...
    handles.sys.mech,...
    handles.idx);


model_ptr = handles.HazOptions.sbh(1);
if ismember(handles.HazOptions.sbh(1),handles.sys.isREG)
    IM_ptr = handles.IM_select.Value;
    handles.site_colors = compute_v(...
        handles.fig,...
        handles.opt,...
        handles.HazOptions,...
        nansum(handles.MRE(:,:,IM_ptr,:,:),4),...
        handles.sys.weight(:,5));
    plot_hazardmap_PSHA(...
        handles.fig,...
        handles.site_colors,...
        handles.HazOptions,...
        handles.h,...
        handles.idx);
end
guidata(hObject,handles)

function addLeg_Callback(hObject, eventdata, handles)

ch = findall(handles.fig,'tag','hazardlegend');
if ~isempty(ch)
    switch ch.Visible
        case 'on',  ch.Visible='off';
        case 'off', ch.Visible='on';
    end
end
guidata(hObject,handles);

function site_menu_psda_Callback(hObject, eventdata, handles)
dsha_lambda2(handles);
guidata(hObject,handles)

function ax2Scale_Callback(hObject, eventdata, handles)
XYSCALE=[handles.ax2.XScale,handles.ax2.YScale];

switch XYSCALE
    case 'linearlinear', handles.ax2.XScale='log';
    case 'loglinear',    handles.ax2.XScale='linear'; handles.ax2.YScale='log';
    case 'linearlog',    handles.ax2.XScale='log';
    case 'loglog',       handles.ax2.XScale='linear'; handles.ax2.YScale='linear';
end

function ax2Grid_Callback(hObject, eventdata, handles)

XY=[handles.ax2.XGrid,handles.ax2.XMinorGrid];
switch XY
    case 'offoff',handles.ax2.XGrid='off'; handles.ax2.XMinorGrid='on';  handles.ax2.YGrid='off'; handles.ax2.YMinorGrid='on';
    case 'offon' ,handles.ax2.XGrid='on';  handles.ax2.XMinorGrid='off'; handles.ax2.YGrid='on';  handles.ax2.YMinorGrid='off';
    case 'onoff' ,handles.ax2.XGrid='on';  handles.ax2.XMinorGrid='on';  handles.ax2.YGrid='on';  handles.ax2.YMinorGrid='on';
    case 'onon'  ,handles.ax2.XGrid='off'; handles.ax2.XMinorGrid='off'; handles.ax2.YGrid='off'; handles.ax2.YMinorGrid='off';
end
guidata(hObject,handles)

function ax2Limits_Callback(hObject, eventdata, handles)
handles.ax2=ax2Control(handles.ax2);

function switchmode_Callback(hObject, eventdata, handles)

switch handles.platformMode
    case 1
        handles.switchmode.CData=handles.form2;
        handles.switchmode.TooltipString='switch to PSHA mode';
        handles.platformMode = 2;
        handles.PSHApannel.Visible='off';
        handles.DSHApannel.Visible='on';
        dsha_lambda2(handles);
        plot_DSHA
        
        
    case 2
        handles.switchmode.CData=handles.form1;
        handles.switchmode.TooltipString='switch to DSHA mode';
        handles.platformMode=1;
        handles.PSHApannel.Visible='on';
        handles.DSHApannel.Visible='off';
        plot_hazard_PSHA(...
            handles.fig,...
            handles.HazOptions,...
            handles.opt,...
            handles.sys.isPCE,...
            handles.MRE,...
            handles.MREPCE,...
            handles.sys.validation,...
            handles.sys.weight(:,5),...
            handles.sys.labelG,...
            handles.sys.mech,...
            handles.idx);
        plot_sites_PSHA
        if ~isempty(handles.h.t)
            handles.po_contours.Enable='on';
            handles.po_contours.Value=1;
            IM_ptr = handles.IM_select.Value;
            handles.site_colors = compute_v(...
                handles.fig,...
                handles.opt,...
                handles.HazOptions,...
                nansum(handles.MRE(:,:,IM_ptr,:,:),4),...
                handles.sys.weight(:,5));
            plot_hazardmap_PSHA(...
                handles.fig,...
                handles.site_colors,...
                handles.HazOptions,...
                handles.h,...
                handles.idx);
        end
        handles.ax2.XScale='log';
        handles.ax2.YScale='log';
end
guidata(hObject,handles)

function SourceLabels_Callback(hObject, eventdata, handles)

val = handles.po_region.Value;
tag = sprintf('sourcetag%g',val);
ch  = findall(handles.ax1,'tag',tag);
switch hObject.Value
    case 0, set(ch,'visible','off')
    case 1, set(ch,'visible','on')
end

function Help_Callback(~,~,~)

function QueryMem_Callback(hObject, eventdata, handles)

M= memory;
flds  = fields(handles);
usedM = zeros(size(flds,1),3);
for i =1:length(flds)
    tmp = handles.(flds{i});
    tmp = whos('tmp');
    sp  = tmp.bytes/M.MemAvailableAllArrays*100;
    if sp>0.01
        usedM(i,:)= [i,tmp.bytes,sp];
    end
end
usedM(usedM(:,1)==0,:)=[];
[~,IND]=sortrows(usedM,-2);
UHand    = M.MemUsedMATLAB/M.MemAvailableAllArrays*100;
tab      = [flds(usedM(IND,1)),num2cell(usedM(IND,2:3))];
Variable = ['Matlab';tab(:,1)];
Size     = [M.MemUsedMATLAB;tab(:,2)];
Usage    = [UHand;tab(:,3)];
T        = table(Variable,Size,Usage);
home
disp(T);

function EnableLightMode_Callback(hObject, eventdata, handles)

switch hObject.Checked
    case 'off'
        hObject.Checked='on';
        handles.engine.Enable='inactive';
    case 'on'
        hObject.Checked='off';
        handles.engine.Enable='on';
end

function About_Callback(~,~,~)
if ispc
    winopen('pshatoolbox_about.txt')
end

function OpenRef_Callback(hObject, eventdata, handles)
if ispc && ~isempty(handles.pdffile) && contains(handles.pdffile,'.pdf')
    winopen(handles.pdffile)
end

function ColorSecondaryLines_Callback(hObject, eventdata, handles)

plot_hazard_PSHA(...
    handles.fig,...
    handles.HazOptions,...
    handles.opt,...
    handles.sys.isPCE,...
    handles.MRE,...
    handles.MREPCE,...
    handles.sys.validation,...
    handles.sys.weight(:,5),...
    handles.sys.labelG,...
    handles.sys.mech,...
    handles.idx);

function restore_axis_Callback(hObject, eventdata, handles)

handles.ax1.XLim=handles.ax1DefaultLimits(1,:);
handles.ax1.YLim=handles.ax1DefaultLimits(2,:);
