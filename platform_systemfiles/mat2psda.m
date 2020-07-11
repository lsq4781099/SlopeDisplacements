function[handles]=mat2psda(handles,varargin)
%#ok<*NODEF>

%% loads sys, opt, h, model
if nargin==3
    pathname=varargin{1};
    filename=varargin{2};
    if contains(filename,'.mat')
        load([pathname,filename],'sys','opt','h')
        handles.sys    = sys; clear sys
        handles.opt    = opt; clear opt
        handles.h      = h;   clear h
    elseif contains(filename,'.txt')
        [handles.sys,handles.opt,handles.h]=loadPSHA(fullfile(pathname,filename));
    end
elseif nargin==5
    handles.sys    = varargin{1};
    handles.opt    = varargin{3};
    handles.h      = varargin{4};
end

[~,kyptr]=intersect(handles.h.param,{'ky','covky'},'stable');
[~,Tsptr]=intersect(handles.h.param,{'Ts','covTs'},'stable');
handles.allky = handles.h.value(:,kyptr);
handles.allTs = handles.h.value(:,Tsptr);

opt.MagDiscrete  = {'uniform',0.1};
n=max(handles.sys.branch(:,3));
for i=1:n
    handles.sys.mrr2(i) = process_truncexp  (handles.sys.mrr2(i) , {'uniform' 0.1});
    handles.sys.mrr3(i) = process_truncnorm (handles.sys.mrr3(i) , {'uniform' 0.1});
    handles.sys.mrr4(i) = process_yc1985    (handles.sys.mrr4(i) , {'uniform' 0.1});
end

%% updates weights
isREG = handles.sys.isREG;
handles.sys.branch         = handles.sys.branch(isREG,:);
handles.sys.branch(:,4)    = handles.sys.branch(:,4)/sum(handles.sys.branch(:,4));

%% loads default models
ME     = handles.ME;
func   = {ME.str};
ptrs   = handles.sys.ptrs(10:14,:);
ptrs   = ptrs-ptrs(1)+1;
[d_default,Sadef,Ddef,SMLIB_default,default_reg]=loadPSDADefaultParam(ME);
txtPSDA = handles.sys.txtPSDA;
% d values
if ~isnan(ptrs(1,1))
    str = txtPSDA(ptrs(1,1):ptrs(1,2),:);
    str = regexp(str,'\ : ','split');
    handles.paramPSDA.d         = eval(str{1}{2});
    handles.paramPSDA.realSa    = str2double(str{2}{2});
    handles.paramPSDA.realD     = str2double(str{3}{2});
    handles.paramPSDA.imhazard  = str{4}{2};
    handles.paramPSDA.rng       = str{5}{2};
    handles.paramPSDA.method    = str{6}{2};
    handles.paramPSDA.optimize  = str{7}{2};
    handles.paramPSDA.kysamples = str2double(str{8}{2});
    handles.paramPSDA.Tssamples = str2double(str{9}{2});
else
    handles.paramPSDA.d         = d_default;
    handles.paramPSDA.realSa    = Sadef;
    handles.paramPSDA.realD     = Ddef;
    handles.paramPSDA.imhazard  = 'full';
    handles.paramPSDA.rng       = 'shuffle';
    handles.paramPSDA.method    = 'MC';
    handles.paramPSDA.optimize  = 'on';
    handles.paramPSDA.kysamples = 0;
    handles.paramPSDA.Tssamples = 0;
end

% displacement model library
if ~isnan(ptrs(2,1))
    str = txtPSDA(ptrs(2,1):ptrs(2,2),:);
    str = regexp(str,'\ handle ','split');
    
    Nmodels = length(str);
    SMLIB(1:Nmodels,1) = struct('id',[],'label',[],'str',[],'func',[],'isregular',[],'integrator',[],'primaryIM',[],'Safactor',[]);
    for i=1:Nmodels
        stri = regexp(str{i}{2},'\ ','split');
        [~,C]=intersect(func,stri{1});
        SMLIB(i).id = str{i}{1};
        SMLIB(i).str        = ME(C).str;
        SMLIB(i).func       = ME(C).func;
        SMLIB(i).isregular  = ME(C).isregular;
        SMLIB(i).integrator = ME(C).integrator;
        SMLIB(i).primaryIM  = ME(C).primaryIM;
        SMLIB(i).Safactor   = ME(C).Safactor;
        SMLIB(i).param      = [];
        if length(stri)>2
            fixparam = struct(stri{2:end});
            flds     = fields(fixparam);
            for j=1:length(flds)
                fixparam.(flds{j})=str2double(fixparam.(flds{j}));
            end
            SMLIB(i).param   = fixparam;
        end
    end
    handles.sys.SMLIB=SMLIB;
else
    handles.sys.SMLIB = SMLIB_default;
end

% Displacement Models for regular PSDA Analysis
if ~isnan(ptrs(3,1))
    str = txtPSDA(ptrs(3,1):ptrs(3,2),:);
    Nmodels = size(str,1)-1;
    T3    = cell(Nmodels,5);
    newline = regexp(str{1},'Weights ','split');
    slopeweights = eval(['[',newline{2},']']);
    slopeweights = slopeweights / sum(slopeweights);
    T3(:,5) = num2cell(slopeweights);
    
    str(1,:)=[];
    
    for j=1:Nmodels
        strj = regexp(str{j},'\ ','split');
        T3{j,1}=strj{1};
        modelassig = struct(strj{2:end});
        T3{j,2} = modelassig.interface;
        T3{j,3} = modelassig.intraslab;
        T3{j,4} = modelassig.crustal;
    end
    handles.T3=T3;
else
    handles.T3=default_reg.T3;
end

%% Displacement Models for PCE PSDA Analysis
if ~isnan(ptrs(4,1))
    str = txtPSDA(ptrs(4,1):ptrs(4,2),:);
    str = regexp(str,'\ ','split');
    Nmodels  = size(str,1);
    DataCDM  = cell(Nmodels,5);
    for j=1:Nmodels
        modelassig   = struct(str{j}{6:end});
        DataCDM{j,1} = str{j}{1};
        DataCDM{j,2} = sprintf('%s,%s,%s',str{j}{3},str{j}{4},str{j}{5});
        DataCDM{j,3} = modelassig.interface;
        DataCDM{j,4} = modelassig.intraslab;
        DataCDM{j,5} = modelassig.crustal;
    end
    isCDMGMM = ~horzcat(SMLIB.isregular);
    handles.tableCDM.ColumnFormat{3}={SMLIB(isCDMGMM).id};
    handles.tableCDM.ColumnFormat{4}={SMLIB(isCDMGMM).id};
    handles.tableCDM.ColumnFormat{5}={SMLIB(isCDMGMM).id};
    handles.tableCDM.Data = DataCDM;
end

%% setup GUI for regular models
if any(isREG)
    handles.pop_site.String=handles.h.id;
    handles.pop_site.Enable='on';
    handles.pop_site.Value=1;
    
    % Table T1
    w1          = handles.sys.branch(:,4);
    id          = compose('Branch %i',isREG');
    handles.T1  = [id,num2cell(w1)];
    
    % Table T2
    kyval    = handles.allky(1,:);
    Tsval    = handles.allTs(1,:);
    handles.T2 = buildPSDA_T2(handles.paramPSDA,kyval,Tsval);
    [handles.tableREG.Data,handles.IJK]=main_psda(handles.T1,handles.T2,handles.T3);
    handles.EditLogicTree.Enable='on';
end


%% ACTIVATE / DEACTIVA REG AND CDM MODES
if ~isempty(handles.tableCDM.Data)
    handles.mode2.Value               = 1;
    handles.mode2.Enable              = 'on';
    handles.tableCDM.Enable           = 'on';
    handles.runCDM.Enable             = 'on';
    handles.CDM_DisplayOptions.Enable = 'on';
else
    handles.mode2.Value               = 0;
    handles.mode2.Enable              = 'off';
    handles.tableCDM.Enable           = 'off';
    handles.runCDM.Enable             = 'inactive';
    handles.CDM_DisplayOptions.Enable = 'inactive';
end

if ~isempty(isREG)
    handles.mode1.Value               = 1;
    handles.mode1.Enable              = 'on';
    handles.tableREG.Enable           = 'on';
    handles.runREG.Enable             = 'on';
    handles.treebutton.Enable         = 'on';
    handles.REG_DisplayOptions.Enable = 'on';
else
    handles.mode1.Value               = 0;
    handles.mode1.Enable              = 'off';
    handles.tableREG.Enable           = 'off';
    handles.runREG.Enable             = 'inactive';
    handles.treebutton.Enable         = 'inactive';
    handles.REG_DisplayOptions.Enable = 'inactive';
end

%% validation data
ind1 =ptrs(5,1);
ind2 =ptrs(5,2);
if ~isnan(ind1)
    line          = regexp(txtPSDA{ind1},'\ ','split');
    handles.sys.D = str2double(line(1,2:end));
    ND            = length(handles.sys.D);
    line          = txtPSDA(ind1+1:ind2,:);
    Nrows         = size(line,1);
    handles.sys.Dlabels     = cell(Nrows,1);
    handles.sys.lambdaDTest = zeros(Nrows,ND);
    
    for i=1:size(line,1)
        line_i = regexp(line{i},'\ ','split');
        lab_i  = strjoin(line_i(1:end-ND),' ');
        line_i = line_i(end-ND+1:end);
        handles.sys.Dlabel{i}=lab_i;
        handles.sys.lambdaDTest(i,:)=str2double(line_i);
    end
end

