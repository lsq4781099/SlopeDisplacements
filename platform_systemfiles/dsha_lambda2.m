function[im,lambdaIS] = dsha_lambda2(handles)

pd = makedist('normal');
sigma = handles.opt.Sigma;
if ~isempty(sigma) && strcmpi(sigma{1},'truncate')
    pd = truncate(pd,-inf,sigma{2});
elseif ~isempty(sigma) && strcmpi(sigma{1},'overwrite')
    pd.sigma=sigma{2};
end

IMptr      = handles.pop_field.Value;
siteptr    = handles.site_menu_psda.Value;
Nsim       = str2double(handles.NumSim.String);
opt        = handles.opt;
im         = handles.opt.im(:,1);
Nim        = size(im,1);
IMs        = [opt.IM1;opt.IM2];
Nsites     = size(handles.h.p,1);
IND        = siteptr + Nsites*(IMptr-1);


[slist,~,sf] = unique(handles.scenarios(:,  1)');
[~,~,Lptr]   = unique(handles.scenarios(:,1:2),'rows','stable');
if size(handles.L,3)==1
    Lptr = ones(size(Lptr));
end
%% EMPIRICAL HAZARD CURVE FROM DSHA ANALYSIS
lambdaIS   = zeros(Nim,1);
Nscen      = size(handles.scenarios,1);
Y_IS       = zeros(Nscen,Nsim);

Nz        = size(handles.L,1);
cont      = 1;

for i=slist
    ind  = find(sf==i)';
    
    for j=ind
        mu  = handles.mulogIM(j,IMptr,siteptr);
        Li  = handles.L(IND,:,Lptr(j));
        Zi  = random(pd,[Nz,Nsim]);
        Y_IS(cont,:) = exp(mu+Li*Zi);
        cont=cont+1;
    end
end
Y_IS    = Y_IS(:);
rate_IS = handles.scenarios(:,9);
rate_IS = repmat(rate_IS/Nsim,1,Nsim);
rate_IS = rate_IS(:);


for i=1:Nim
    lambdaIS(i)=rate_IS'*(Y_IS>im(i));
end

%% kmean Clusters
if ~isempty(handles.krate)
    lambdaKM = zeros(Nim,1);
    kY       = handles.kY(:,IND);
    krate    = handles.krate;
    for i=1:Nim
        lambdaKM(i)=sum(krate(kY>=im(i)));
    end
end

%% SEISMIC HAZARD CURVE FROM PSHA ANALYSIS (FOR COMPARISON PURPOSES)
IM         = IMs(IMptr);
h.p        = handles.h.p(siteptr,:);
h.param    = handles.h.param;
h.value    = handles.h.value(siteptr,:);

Nsource    = length(handles.source);
MRE        = runhazard1(im,IM,h,opt,handles.source,Nsource,1);
lambdaPSHA = permute(MRE,[2,4,1,3]);
lambdaPSHA = nansum(lambdaPSHA,2);


%% plotting
delete(findall(handles.ax2,'Type','line'));
handles.ax2.ColorOrderIndex=1;
if ~isempty(handles.krate)
    plot(handles.ax2,im,lambdaPSHA ,'k--' ,'tag','lambda_empirical','ButtonDownFcn',{@myfun,handles}),
    plot(handles.ax2,im,lambdaIS   ,'.-' ,'tag','lambda_empirical','ButtonDownFcn',{@myfun,handles}),
    plot(handles.ax2,im,lambdaKM   ,'.-' ,'tag','lambda_empirical','ButtonDownFcn',{@myfun,handles}),
    LEG={'PSHA',sprintf('DSHA (IS %g)',length(rate_IS)),sprintf('DSHA (kM %g)',length(krate))};
else
    plot(handles.ax2,im,lambdaPSHA ,'k--' ,'tag','lambda_empirical','ButtonDownFcn',{@myfun,handles}),    
    plot(handles.ax2,im,lambdaIS   ,'.-' ,'tag','lambda_empirical','ButtonDownFcn',{@myfun,handles}),
    LEG={'PSHA',sprintf('DSHA (IS %g)',length(rate_IS))};
end
L=legend(handles.ax2,LEG);
L.Visible='on';
L.Location='northeast';
L.Box='off';



%% puts data in clipboard
if isempty(handles.krate)
    data  = num2cell([zeros(1,3);[im,lambdaIS,lambdaIS]]);
    data{1,1}='im';
    data{1,2}='PSHA';
    data{1,3}='IS';
else
    data  = num2cell([zeros(1,4);[im,lambdaIS,lambdaIS,lambdaKM]]);
    data{1,1}='im';
    data{1,2}='PSHA';
    data{1,3}='IS';
    data{1,4}='KM';
end
c2 = uicontextmenu;
uimenu(c2,'Label','Copy data','Callback'            ,{@data2clipboard_uimenu,data(2:end,:)});
uimenu(c2,'Label','Copy data & headers','Callback'  ,{@data2clipboard_uimenu,data(1:end,:)});
uimenu(c2,'Label','Undock','Callback'               ,{@figure2clipboard_uimenu,handles.ax2});
set(handles.ax2,'uicontextmenu',c2);




function[]=myfun(hObject, eventdata, handles) %#ok<INUSL>

H=datacursormode(handles.fig);
set(H,'enable','on','DisplayStyle','window','UpdateFcn',{@gethazarddata,handles.HazOptions.dbt(1)});
w = findobj('Tag','figpanel');
set(w,'Position',[ 409   485   150    60]);

function output_txt = gethazarddata(~,event_obj,dbt)

pos = get(event_obj,'Position');

if dbt==1
    output_txt = {...
        ['IM   : ',num2str(pos(1),4)],...
        ['Rate : ',num2str(pos(2),4)],...
        ['T    : ',num2str(1/pos(2),4),' years']};
end

if dbt==0
    output_txt = {...
        ['IM             : ',num2str(pos(1),4)],...
        ['P(IM>im|t) : ',num2str(pos(2),4)]};
end
