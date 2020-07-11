function [handles]=GMMValidation_CY2008(handles,filename)

switch filename
    
    case 'CY2008_1.png'
        %M, Rrup,    Rjb,   Rx, Ztor,    dip, lambda,   event, Vs30, Z10, FVS30
        handles.e1.String=6.8;    % M
        handles.e2.String=5;      % Rrup
        handles.e3.String=0;      % Rjb
        handles.e4.String=5;      % Rx
        handles.e5.String=0;      % Ztor
        handles.e6.String=60;     % dip
        handles.e7.String=1100;   % Vs30
        handles.e8.String=5.744229155; %Z1.0
        handles.e9.Value=4;         % mechanism
        handles.e10.Value=1;         % event
        handles.e11.Value=2;        %FVS30
        plotgmpe(handles);
        handles.ax1.XLim=[0.01  10];
        handles.ax1.YLim=[0.001 10];
        handles.ax1.XScale='log';
        handles.ax1.YScale='log';
        
    case 'CY2008_2.png'
        %M, Rrup,    Rjb,   Rx, Ztor,    dip, lambda,   event, Vs30, Z10, FVS30
        handles.e2.String=1;      % Rrup
        handles.e3.String=0;      % Rjb
        handles.e4.String=1;      % Rx
        handles.e5.String=0;      % Ztor
        handles.e6.String=90;     % dip
        handles.e7.String=520;    % Vs30
        handles.e8.String=exp(28.5-3.82/8*log(520^8+378.8^8)); %Z1.0
        handles.e9.Value=1;         % mechanism
        handles.e10.Value=1;         % event
        handles.e11.Value=2;           %FVS30
        handles.e1.String=5.5;   plotgmpe(handles);
        handles.e1.String=6.5;   plotgmpe(handles);
        handles.e1.String=7.5;   plotgmpe(handles);
        handles.e1.String=8.5;   plotgmpe(handles);
        handles.ax1.XLim=[0.01  10];
        handles.ax1.YLim=[0.0007 2];
        handles.ax1.XScale='log';
        handles.ax1.YScale='log';
        
    case 'CY2008_3.png'
        %M, Rrup,    Rjb,   Rx, Ztor,    dip, lambda,   event, Vs30, Z10, FVS30
        handles.e2.String=10;     % Rrup
        handles.e3.String=0;      % Rjb
        handles.e4.String=10;     % Rx
        handles.e5.String=0;      % Ztor
        handles.e6.String=90;     % dip
        handles.e7.String=520;    % Vs30
        handles.e8.String=exp(28.5-3.82/8*log(520^8+378.8^8)); %Z1.0
        handles.e9.Value=1;       % mechanism
        handles.e10.Value=1;      % event
        handles.e11.Value=2;      %FVS30
        handles.e1.String=5.5;   plotgmpe(handles);
        handles.e1.String=6.5;   plotgmpe(handles);
        handles.e1.String=7.5;   plotgmpe(handles);
        handles.e1.String=8.5;   plotgmpe(handles);
        handles.ax1.XLim=[0.01  10];
        handles.ax1.YLim=[0.0007 2];
        handles.ax1.XScale='log';
        handles.ax1.YScale='log';
        
    case 'CY2008_4.png'
        %M, Rrup,    Rjb,   Rx, Ztor,    dip, lambda,   event, Vs30, Z10, FVS30
        handles.e2.String=30;     % Rrup
        handles.e3.String=0;      % Rjb
        handles.e4.String=30;     % Rx
        handles.e5.String=0;      % Ztor
        handles.e6.String=90;     % dip
        handles.e7.String=520;    % Vs30
        handles.e8.String=exp(28.5-3.82/8*log(520^8+378.8^8)); %Z1.0
        handles.e9.Value=1;       % mechanism
        handles.e10.Value=1;      % event
        handles.e11.Value=2;      %FVS30
        handles.e1.String=5.5;   plotgmpe(handles);
        handles.e1.String=6.5;   plotgmpe(handles);
        handles.e1.String=7.5;   plotgmpe(handles);
        handles.e1.String=8.5;   plotgmpe(handles);
        handles.ax1.XLim=[0.01  10];
        handles.ax1.YLim=[0.0007 2];
        handles.ax1.XScale='log';
        handles.ax1.YScale='log';
        
    case 'CY2008_5.png'
        %M, Rrup,    Rjb,   Rx, Ztor,    dip, lambda,   event, Vs30, Z10, FVS30
        handles.e2.String=70;      % Rrup
        handles.e3.String=0;      % Rjb
        handles.e4.String=70;      % Rx
        handles.e5.String=0;      % Ztor
        handles.e6.String=90;     % dip
        handles.e7.String=520;    % Vs30
        handles.e8.String=num2str(exp(28.5-3.82/8*log(520^8+378.8^8))); %Z1.0
        handles.e9.Value=1;         % mechanism
        handles.e10.Value=1;         % event
        handles.e11.Value=2;           %FVS30
        handles.e1.String=5.5;   plotgmpe(handles);
        handles.e1.String=6.5;   plotgmpe(handles);
        handles.e1.String=7.5;   plotgmpe(handles);
        handles.e1.String=8.5;   plotgmpe(handles);
        handles.ax1.XLim=[0.01  10];
        handles.ax1.YLim=[0.0007 2];
        handles.ax1.XScale='log';
        handles.ax1.YScale='log';
        
end