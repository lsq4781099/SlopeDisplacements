function plotgmpe(handles,amp)

if handles.HoldPlot.Value==0
    delete(findall(handles.ax1,'type','line'))
    handles.AutoScale.Value=1;
    handles.ax1.XLimMode='auto';
    handles.ax1.YLimMode='auto';
    handles.ax1.XTickMode='auto';
    handles.ax1.YTickMode='auto';
    handles.ax1.XTickLabelMode='auto';
    handles.ax1.YTickLabelMode='auto';
end

epsilon = handles.epsilon;
Neps    = length(epsilon);

% gmpe parameters
param     = mGMPEgetparam(handles);

if nargin==1
    amp=1;
end

if handles.rad1.Value==1
    IM=real(handles.IM); 
    LB=imag(handles.IM);
    
    LB=LB(IM>=0);
    IM=IM(IM>=0);
    
    uLB = unique(LB);
   
    if length(uLB)>1
        str0 = {'SA(T)','SV(T)','SD(T)','H/V(T)'};
        str0 = str0(uLB+1);
        [indx,tf] = listdlg('ListString',str0,'SelectionMode','single');
        if tf==0
            uLB=uLB(1);
        else
            uLB=uLB(indx);
        end
    end
    
    switch uLB
        case 0, YLABEL = 'SA(g)';    IM = IM(LB==0); LB = LB(LB==0);
        case 1, YLABEL = 'SV(cm/s)'; IM = IM(LB==1); LB = LB(LB==1);
        case 2, YLABEL = 'SD(cm)';   IM = IM(LB==2); LB = LB(LB==2);
        case 3, YLABEL = 'H/V';      IM = IM(LB==3); LB = LB(LB==3);
    end
    
    Sa        = nan(Neps,length(IM)+1);
    for i=1:length(IM)
        To = IM(i)+LB(i)*1i;
        [lny,sigma] = handles.fun(To,param{:});
        for j=1:Neps
            Sa(j,i) = exp(lny+epsilon(j)*sigma);
        end
    end
    
    switch handles.HoldPlot.Value
        case 0
            ch=findall(handles.ax1,'tag','curves');delete(ch);
            handles.ax1.ColorOrderIndex = 1;
        case 1
    end
    x = repmat([IM,nan],1,Neps);
    
    Sa = Sa';
    y = Sa(:)';
    plot(handles.ax1,x,amp*y,'tag','curves','ButtonDownFcn',@click_on_curve,'linewidth',1);
    handles.xlabel=xlabel(handles.ax1,'T(s)','fontsize',10);
    handles.ylabel=ylabel(handles.ax1,YLABEL,'fontsize',10);
else
    
    imptr  = handles.targetIM.Value;
    [Rrup,param2] = mGMPErrupLoop(handles.fun,param,handles.SUB,handles.SC);
    if isempty(param2)
        Sa = nan(length(Rrup)+1,Neps);
    else
        IM     = handles.IM(imptr);
        [lny,sigma] = handles.fun(IM,param2{:});
        Sa = nan(length(lny)+1,Neps);
        for i=1:Neps
            Sa(1:end-1,i) = exp(lny+epsilon(i)*sigma);
        end
    end
    
    switch handles.HoldPlot.Value
        case 0
            ch=findall(handles.ax1,'tag','curves');delete(ch);
            handles.ax1.ColorOrderIndex = 1;
        case 1
    end
    x = repmat([Rrup;nan],Neps,1);
    y = Sa(:);
    plot(handles.ax1,x,amp*y,'tag','curves','ButtonDownFcn',@click_on_curve,'linewidth',1);
    xlabel(handles.ax1,'Rrup(km)','fontsize',10);
    ylabel(handles.ax1,addIMunits(IM),'fontsize',10);
end

% updates legend
ch=findall(handles.figure1,'type','legend');
new_str = handles.GMPEselect.String{handles.GMPEselect.Value};

switch handles.HoldPlot.Value
    case 0
        delete(ch);
        ch=legend(handles.ax1,new_str);
        ch.Box='off';
    case 1
        ch.String{end}=new_str;
end

switch handles.disp_legend.Value
    case 0, ch.Visible='off';
    case 1, ch.Visible='on';
end



function click_on_curve(hObject,eventdata)

switch eventdata.Button
    case 1 %click izquierdo
    case 3 %click derecho
        delete(hObject)
end

