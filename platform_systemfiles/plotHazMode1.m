function str=plotHazMode1(fig,haz,opt,isPCE,MRE,MREPCE,validation,weights,idx)

ax1          = findall(fig,'tag','ax1');
ax2          = findall(fig,'tag','ax2');
OpenRef      = findall(fig,'tag','OpenRef');
site_menu    = findall(fig,'tag','site_menu'); 
IM_select    = findall(fig,'tag','IM_select');

if isempty(idx)
    site       = site_menu.Value;
else
    site       = idx(site_menu.Value);
end
IM_ptr     = IM_select.Value;

%% -------------  initialization ----------------------------------
if size(opt.im,2)==1
    im  = opt.im';
else
    im  = opt.im(:,IM_ptr)';
end

delete(findall(ax1,'tag','scenario'));
delete(findall(ax2,'type','line'));
delete(findall(ax2,'tag','histogram'));
ax2.ColorOrderIndex = 1;

%% -------------  compute lambda1----------------------------------
lambda1  = nansum(MRE(site,:,IM_ptr,:,:),4);
lambda1  = permute(lambda1,[5,2,1,3,4]);

for kk=isPCE
    lambda1PCE = nansum(MREPCE{kk}(site,:,IM_ptr,:,:),4);
    switch haz.pce(1)
        case 1, lambda1(kk,:) = real(exp(nanmean(log(lambda1PCE),5))); % trucazo
        case 0, lambda1(kk,:) = prctile(lambda1PCE,haz.pce(3),5);
    end
end

%% -------------  compute lambda0----------------------------------
switch find(haz.avg(1:3))
    case 1
        lambda0 = prod(bsxfun(@power,lambda1,weights(:)),1);
    case 2
        weights  = haz.rnd;
        lambda0 = prod(bsxfun(@power,lambda1,weights(:)),1);
    case 3
        lambda0 = prctile(lambda1,haz.avg(4),1);
end

if haz.dbt(1)==1
    y0 = lambda0;
    y1 = lambda1;
end

if haz.dbt(3)==1
    y0 = 1-exp(-lambda0*haz.dbt(4));
    y1 = 1-exp(-lambda1*haz.dbt(4));
end

% notnan = true(size(im));
if ~isempty(validation)
    lambdaTest=validation(site+1,:);
    IM       = validation(1,:);
    switch haz.dbt(3)
        case 0
            yT = lambdaTest;
        case 1
            yT = 1-exp(-lambdaTest*haz.dbt(4));
    end
end

%% -------------  plot hazard--------------------------------------
y0clean = y0; y0clean(y0<0)= nan;
y1clean = y1; y1clean(y1<0)=nan;

switch haz.avg(5)
    case 0
        plot(ax2,im',y0clean','.-','ButtonDownFcn',{@myfun,fig,haz},'tag','lambda0');
        if ~isempty(validation)
            plot(ax2,IM',yT','ko','tag','lambdaTest','markerfacecolor','none');
        end
    case 1
        ax2.ColorOrderIndex=2;
        plot(ax2,im',y1clean','-','ButtonDownFcn',{@myfun,fig,haz},'tag','lambda1','visible','on')
        plot(ax2,im',y0clean','.-','color',[0 0.447 0.741],'ButtonDownFcn',{@myfun,fig,haz},'tag','lambda0');
        if ~isempty(validation)
            plot(ax2,IM',yT','ko','tag','lambdaTest','markerfacecolor','none');
        end
end

%% -------------  xlabel & ylabel ----------------------------------
if iscell(IM_select.String)
    xlabel(ax2,addIMunits(IM_select.String{IM_ptr}),'fontsize',8)
else
    xlabel(ax2,addIMunits(IM_select.String),'fontsize',8)
end

switch find(haz.dbt(1:3))
    case 1,ylabel(ax2,'Mean Rate of Exceedance','fontsize',8)
    case 2,ylabel(ax2,'Mean Rate Density','fontsize',8)
    case 3,ylabel(ax2,'Probability of Exceedance','fontsize',8)
end

%% -------------  legend  -----------------------------------------
switch find(haz.avg(1:3))
    case 1, str1 = 'Default Weights';
    case 2, str1 = 'Random Weights';
    case 3, str1 = sprintf('Percentile %g',haz.avg(4));
end
if haz.avg(5)
    str2 = compose('Branch %i',1:size(MRE,5));
else
    str2={};
end
if ~isempty(validation)
    str3 = {'Benchmark'};
else
    str3 = {};
end
str = [str2,str1,str3];

% -------------  ui context menu ------------------------------------------
IMstr = IM_select.String{IM_ptr};
notnan=find(~isnan(sum(lambda1,2)));

if ~haz.avg(5)
    data  = num2cell([zeros(1,2);[im',y0']]);
else
    Ncol  = size(y0,1)+1+length(notnan);
    data  = num2cell([zeros(1,Ncol);[im',y1(notnan,:)',y0']]);
end


data{1,1}=IMstr;
if ~isempty(validation)
    if size(data,2)==2
        data(1,2:end)={str};
    else
        data(1,2:end)=str;
    end
    OpenRef.Visible='off';
else
    data(1,2:end)=str;
    OpenRef.Visible='on';
end

c2 = uicontextmenu;
uimenu(c2,'Label','Copy data','Callback'            ,{@data2clipboard_uimenu,data(2:end,:)});
uimenu(c2,'Label','Copy data & headers','Callback'  ,{@data2clipboard_uimenu,data(1:end,:)});
uimenu(c2,'Label','Undock','Callback'               ,{@figure2clipboard_uimenu,ax2});
uimenu(c2,'Label','Undock & compare','Callback'     ,{@figurecompare_uimenu,ax2});
set(ax2,'uicontextmenu',c2);

function[]=myfun(hObject, eventdata, fig,haz) %#ok<INUSL>

H=datacursormode(fig);
set(H,'enable','on','DisplayStyle','window','UpdateFcn',{@gethazarddata,haz.dbt(1)});
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
