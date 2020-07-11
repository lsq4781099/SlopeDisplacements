function str=plotHazMode3(fig,haz,opt,MREPCE,idx)

ax1          = findall(fig,'tag','ax1');
ax2          = findall(fig,'tag','ax2');
site_menu    = findall(fig,'tag','site_menu'); 
IM_select    = findall(fig,'tag','IM_select');

site   = site_menu.Value;
IM_ptr = IM_select.Value;

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

%% -------------  compute lambda1 and lambda0 ----------------------------
model_ptr = haz.sbh(1);
lambda1   = MREPCE{model_ptr}(site,:,IM_ptr,:,:);
lambda1   = nansum(lambda1,4);
lambda1   = permute(lambda1,[5,2,1,3,4]);

if haz.dbt(2)==1
    Nreal   = size(lambda1,1);
    for kkk=1:Nreal
        lambda1(kkk,:)= MREDer(im,-lambda1(kkk,:),1)';
    end
    lambda1(lambda1<0)=0;
end

switch haz.pce(1)
    case 1
        loglambda1=log(lambda1);
        loglambda1(isinf(loglambda1))=nan;
        lambda0 = real(exp(nanmean(loglambda1,1))); % trucazo
    case 0 
        lambda0 = prctile(lambda1,haz.pce(3),1);
end

switch haz.dbt(3)
    case 0
        y0 = lambda0;
        y1 = lambda1;
    case 1
        y0 = 1-exp(-lambda0*haz.dbt(4));
        y1 = 1-exp(-lambda1*haz.dbt(4));
end

%% -------------  plot hazard--------------------------------------
y0(y0<0)= nan;
y1(y1<0)=nan;


if haz.sbh(4)
    plot(ax2,nan,nan,'-','color',[1 1 1]*0.8,'ButtonDownFcn',{@myfun,fig,haz},'tag','lambda1','visible','on')
    plot(ax2,im',y1','-','color',[1 1 1]*0.8,'ButtonDownFcn',{@myfun,fig,haz},'tag','lambda1','visible','on','HandleVisibility','off')
end
ax2.ColorOrderIndex=1;
plot(ax2,im',y0','-','ButtonDownFcn',{@myfun,fig,haz},'tag','lambda0','linewidth',2);
ax2.Layer='top';

%% -------------  xlabel & ylabel ----------------------------------
if iscell(IM_select.String)
    xlabel(ax2,addIMunits(IM_select.String{IM_ptr}),'fontsize',8)
else
    xlabel(ax2,addIMunits(IM_select.String),'fontsize',8)
end

switch find(haz.dbt(1:3))
    case 1,ylabel(ax2,'Mean Rate of Exceedance','fontsize',10)
    case 2,ylabel(ax2,'Mean Rate Density','fontsize',10)
    case 3,ylabel(ax2,'Probability of Exceedance','fontsize',10)
end

%% -------------  legend  -----------------------------------------
switch opt.PCE{2}
    case 'PC'
        str = sprintf('Branch %g (PCE)',model_ptr);
        if haz.sbh(4), str = [{'Simulations (PCE)'},str];  end
    case 'MC'
        str = sprintf('Branch %g (MCS)',model_ptr);
        if haz.sbh(4), str = [{'Simulations (MCS)'},str];  end
end

%%  -------------  ui context menu ------------------------------------------
IMstr = IM_select.String{IM_ptr};
data  = num2cell([zeros(1,2);[im',y0']]);

data{1,1}=IMstr;
if size(data,2)==2
    data(1,2:end)={str};
else
    data(1,2:end)=str;
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
