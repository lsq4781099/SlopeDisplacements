function [handles]=GMMValidation_BA2008(handles,filename)

switch filename
    case 'BA2008_1.png'
        %M, Rjb, Fault_Type, Vs30
        handles.e2.String=0;    % Rjb
        handles.e3.String=760;  % Vs30
        handles.e4.Value=1;     % Fault_type
        handles.e1.String=8; plotgmpe(handles,980.66);
        handles.e1.String=7; plotgmpe(handles,980.66);
        handles.e1.String=6; plotgmpe(handles,980.66);
        handles.e1.String=5; plotgmpe(handles,980.66);
        
        handles.ax1.XLim=[0.01  10];
        handles.ax1.YLim=[0.01 2000];
        handles.ax1.XScale='log';
        handles.ax1.YScale='log';
        
    case 'BA2008_2.png'
        %M, Rjb, Fault_Type, Vs30
        handles.e2.String=30;   % Rjb
        handles.e3.String=760;  % Vs30
        handles.e4.Value=1;     % Fault_type
        handles.e1.String=8; plotgmpe(handles,980.66);
        handles.e1.String=7; plotgmpe(handles,980.66);
        handles.e1.String=6; plotgmpe(handles,980.66);
        handles.e1.String=5; plotgmpe(handles,980.66);
        
        handles.ax1.XLim=[0.01  10];
        handles.ax1.YLim=[0.01 2000];
        handles.ax1.XScale='log';
        handles.ax1.YScale='log';
        
    case 'BA2008_3.png'
        %M, Rjb, Fault_Type, Vs30
        handles.e2.String=200;   % Rjb
        handles.e3.String=760;   % Vs30
        handles.e4.Value=1;      % Fault_type
        handles.e1.String=8; plotgmpe(handles,980.66);
        handles.e1.String=7; plotgmpe(handles,980.66);
        handles.e1.String=6; plotgmpe(handles,980.66);
        handles.e1.String=5; plotgmpe(handles,980.66);
        
        handles.ax1.XLim=[0.01  10];
        handles.ax1.YLim=[0.01 2000];
        handles.ax1.XScale='log';
        handles.ax1.YScale='log';
end