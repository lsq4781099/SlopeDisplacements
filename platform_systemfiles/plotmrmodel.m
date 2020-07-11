function plotmrmodel(handles)


if handles.MRselect.Value==5
    NMmin=regexp(handles.e3.String,'\ ','split');
    NMmin=str2double(NMmin{1});
else
    NMmin=str2double(handles.e1.String);
end

% gmpe parameters
[~,param] = mMRgetparam(handles);
x         = linspace(5,10,200);
TOL       = 1e-10;
switch func2str(handles.fun)
    case 'truncexp'
        Mmin  = param(2);
        Mmax  = param(3);
        tr1=-TOL + Mmin;
        tr2=-TOL + Mmax;
        tr3=+TOL + Mmax;
        x = sort([x,tr1,Mmin,Mmax,tr2,tr3]);
        x(x>Mmax)=NaN;
    case 'truncnorm'
        Mchar   = param(3);
        sigmaM  = param(4);
        tr1=-TOL + Mchar-sigmaM;
        tr2=       Mchar-sigmaM;
        tr3=+TOL + Mchar+sigmaM;
        tr4=       Mchar+sigmaM;
        x = sort([x,tr1,tr2,tr3,tr4]);
        x(x>param(2))=NaN;
    case 'delta'
        M   = param(1);
        tr1=-TOL + M;
        tr2=+TOL + M;
        x = sort([x,tr1,M,tr2]);
        
    case 'yc1985'
        Mmin   = param(2);
        Mchar  = param(3);
        tr1 = -TOL + Mmin;
        tr2 = Mmin;
        tr3 = Mchar-0.25-TOL;
        tr4 = Mchar-0.25;
        tr5 = Mchar-0.25+TOL;
        tr6 = Mchar+0.25;
        tr7 = Mchar+0.25+TOL;
        x = sort([x,tr1,tr2,tr3,tr4,tr5,tr6,tr7]);
        x(x>Mchar+0.25+TOL)=NaN;
    case 'magtable'
        Mmin       = param(1);
        binwidth   = param(2);
        NMmin      = param(3);
        occurrance = param(3:end);
        nM         = length(occurrance);
        x          = Mmin:binwidth:(Mmin+binwidth*(nM-1));
end

switch handles.uibuttongroup1.SelectedObject.String
    case 'rate',[~,y] = handles.fun(x,param); ylabel(handles.ax1,'Magnitude rate');y = NMmin*(1-y);
    case 'pdf', [y,~] = handles.fun(x,param); ylabel(handles.ax1,'Magnitude pdf');
    case 'cdf', [~,y] = handles.fun(x,param); ylabel(handles.ax1,'Magnitude cdf');
end

switch handles.HoldPlot.Value
    case 0
        ch=findall(handles.ax1,'tag','curves');delete(ch);
        handles.ax1.ColorOrderIndex = 1;
    case 1
end
y(y<0)=nan;
plot(handles.ax1,x,y,'-','tag','curves','ButtonDownFcn',@click_on_curve,'handlevisibility','off');
handles.xlabel=xlabel(handles.ax1,'Magnitude','fontsize',10);


%% -------------- data2clipboard ---------------------------
% cF=get(0,'format');
% format long g
% data  = num2cell([x(:),y(:)]);
% c2 = uicontextmenu;
% uimenu(c2,'Label','Copy data','Callback'    ,{@data2clipboard_uimenu,data});
% uimenu(c2,'Label','Undock','Callback'       ,{@figure2clipboard_uimenu,handles.ax1});
% set(handles.ax1,'uicontextmenu',c2);
% format(cF);


function click_on_curve(hObject,eventdata)

switch eventdata.Button
    case 1 %click izquierdo
    case 3 %click derecho
        delete(hObject)
end

