function [handles]=GMMValidation_AS1997h(handles,filename)

e1  = handles.e1;
e2  = handles.e2;
e3  = handles.e3;
e4  = handles.e4;
e5  = handles.e5;
e6  = handles.e6;


switch filename
    case 'AS1997h_1.png'
        e1.String = 7;
        e4.Value  = 3; % {'reverse','reverse-oblique','other'};
        e5.Value  = 2; % {'hangingwall','footwall'};
        e6.Value  = 1; % {'measured','inferred'};
        e2.String = 1;   e3.Value = 1; plotgmpe(handles);
        e2.String = 10;  e3.Value = 1; plotgmpe(handles);
        e2.String = 30;  e3.Value = 1; plotgmpe(handles);
        e2.String = 100; e3.Value = 1; plotgmpe(handles);
        
        e2.String = 1;   e3.Value = 2; plotgmpe(handles);
        e2.String = 10;  e3.Value = 2; plotgmpe(handles);
        e2.String = 30;  e3.Value = 2; plotgmpe(handles);
        e2.String = 100; e3.Value = 2; plotgmpe(handles);
        
        
        handles.ax1.XLim=[0.01  10];
        handles.ax1.YLim=[0.01 10];
        handles.ax1.XScale='log';
        handles.ax1.YScale='log';
        
    case 'AS1997h_2.png'
        e2.String = 10;
        e4.Value  = 3; % {'reverse','reverse-oblique','other'};
        e5.Value  = 2; %{'hangingwall','footwall'};
        e6.Value  = 1; % {'measured','inferred'};
        
        e1.String = 5.5;  e3.Value = 1; plotgmpe(handles);
        e1.String = 6.5;  e3.Value = 1; plotgmpe(handles);
        e1.String = 7.5;  e3.Value = 1; plotgmpe(handles);
        
        e1.String = 5.5;  e3.Value = 2; plotgmpe(handles);
        e1.String = 6.5;  e3.Value = 2; plotgmpe(handles);
        e1.String = 7.5;  e3.Value = 2; plotgmpe(handles);
        
        handles.ax1.XLim=[0.01  10];
        handles.ax1.YLim=[0.01 10];
        handles.ax1.XScale='log';
        handles.ax1.YScale='log'; 
end

handles.e1  = e1;
handles.e2  = e2;
handles.e3  = e3;
handles.e4  = e4;
handles.e5  = e5;
handles.e6  = e6;
