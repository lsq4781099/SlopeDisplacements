function [handles]=GMMValidation_CB2008(handles,filename)


switch filename
    
    case 'CB2008_1.png'
        %M, Rrup, Rjb, Ztor, dip, mechanism, Vs30, Zvs, arb, A1100
        handles.e2.String=0;      % Rrup
        handles.e3.String=0;      % Rjb
        handles.e4.String=0;      % Ztor
        handles.e5.String=90;     % dip
        handles.e6.String=760;    % Vs30
        handles.e7.String=2;      % 2
        handles.e8.Value=1;         % mechanism
        handles.e9.Value=1;       %FVS30
        
        handles.e1.String=5;   plotgmpe(handles);
        handles.e1.String=6;   plotgmpe(handles);
        handles.e1.String=7;   plotgmpe(handles);
        handles.e1.String=8;   plotgmpe(handles);
        handles.ax1.XLim=[0.01  10];
        handles.ax1.YLim=[1e-5 1.001];
        handles.ax1.XScale='log';
        handles.ax1.YScale='log';
        
    case 'CB2008_2.png'
        %M, Rrup, Rjb, Ztor, dip, mechanism, Vs30, Zvs, arb, A1100
        handles.e2.String=10;     % Rrup
        handles.e3.String=0;      % Rjb
        handles.e4.String=0;      % Ztor
        handles.e5.String=90;     % dip
        handles.e6.String=760;    % Vs30
        handles.e7.String=2;      % 2
        handles.e8.Value=1;       % mechanism
        handles.e9.Value=1;       %FVS30
        
        handles.e1.String=5;   plotgmpe(handles);
        handles.e1.String=6;   plotgmpe(handles);
        handles.e1.String=7;   plotgmpe(handles);
        handles.e1.String=8;   plotgmpe(handles);
        handles.ax1.XLim=[0.01  10];
        handles.ax1.YLim=[1e-5 1.001];
        handles.ax1.XScale='log';
        handles.ax1.YScale='log';
        
    case 'CB2008_3.png'
        %M, Rrup, Rjb, Ztor, dip, mechanism, Vs30, Zvs, arb, A1100
        handles.e2.String=50;     % Rrup
        handles.e3.String=0;      % Rjb
        handles.e4.String=0;      % Ztor
        handles.e5.String=90;     % dip
        handles.e6.String=760;    % Vs30
        handles.e7.String=2;      % 2
        handles.e8.Value=1;       % mechanism
        handles.e9.Value=1;       %FVS30
        handles.e10.String=0;     %A1100
        handles.e1.String=5;   plotgmpe(handles);
        handles.e1.String=6;   plotgmpe(handles);
        handles.e1.String=7;   plotgmpe(handles);
        handles.e1.String=8;   plotgmpe(handles);
        handles.ax1.XLim=[0.01  10];
        handles.ax1.YLim=[1e-5 1.001];
        handles.ax1.XScale='log';
        handles.ax1.YScale='log';
        
    case 'CB2008_4.png'
        %M, Rrup, Rjb, Ztor, dip, mechanism, Vs30, Zvs, arb, A1100
        handles.e2.String=200;    % Rrup
        handles.e3.String=0;      % Rjb
        handles.e4.String=0;      % Ztor
        handles.e5.String=90;     % dip
        handles.e6.String=760;    % Vs30
        handles.e7.String=2;      % 2
        handles.e8.Value=1;       % mechanism
        handles.e9.Value=1;       % FVS30
        
        handles.e1.String=5;   plotgmpe(handles);
        handles.e1.String=6;   plotgmpe(handles);
        handles.e1.String=7;   plotgmpe(handles);
        handles.e1.String=8;   plotgmpe(handles);
        handles.ax1.XLim=[0.01  10];
        handles.ax1.YLim=[1e-5 1.001];
        handles.ax1.XScale='log';
        handles.ax1.YScale='log';
        
end