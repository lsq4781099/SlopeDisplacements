function [handles]=GMMValidation_BCHydro2012(handles,filename)

switch filename
    case 'BCHydro2012_1.png'
        % Figure 10a from Abrahamson 2016
        handles.e2.String=25;
        handles.e3.String=[];
        handles.e4.String=760;
        handles.e5.Value=1;
        handles.e6.Value=1;
        handles.e7.Value=2;
        
        handles.e1.String=7.0; plotgmpe(handles);
        handles.e1.String=8.0; plotgmpe(handles);
        handles.e1.String=9.0; plotgmpe(handles);
        handles.ax1.XLim=[0.01 10];
        handles.ax1.YLim=[0.001 10];
        handles.ax1.XScale='log';
        handles.ax1.YScale='log';
        
    case 'BCHydro2012_2.png'
        % Figure 10b from Abrahamson 2016
        handles.e2.String=50;
        handles.e3.String=[];
        handles.e4.String=760;
        handles.e5.Value=1;
        handles.e6.Value=1;
        handles.e7.Value=2;
        
        handles.e1.String=7.0; plotgmpe(handles);
        handles.e1.String=8.0; plotgmpe(handles);
        handles.e1.String=9.0; plotgmpe(handles);
        handles.ax1.XLim=[0.01 10];
        handles.ax1.YLim=[0.001 10];
        handles.ax1.XScale='log';
        handles.ax1.YScale='log';
        
    case 'BCHydro2012_3.png'
        % Figure 10c from Abrahamson 2016
        handles.e2.String=100;
        handles.e3.String=[];
        handles.e4.String=760;
        handles.e5.Value=1;
        handles.e6.Value=1;
        handles.e7.Value=2;        
        
        handles.e1.String=7.0; plotgmpe(handles);
        handles.e1.String=8.0; plotgmpe(handles);
        handles.e1.String=9.0; plotgmpe(handles);
        handles.ax1.XLim=[0.01 10];
        handles.ax1.YLim=[0.0001 1];
        handles.ax1.XScale='log';
        handles.ax1.YScale='log';
        
    case 'BCHydro2012_4.png'
        % Figure 10d from Abrahamson 2016
        handles.e2.String=200;
        handles.e3.String=[];
        handles.e4.String=760;
        handles.e5.Value=1;
        handles.e6.Value=1;
        handles.e7.Value=2;   
        
        handles.e1.String=7.0; plotgmpe(handles);
        handles.e1.String=8.0; plotgmpe(handles);
        handles.e1.String=9.0; plotgmpe(handles);
        handles.ax1.XLim=[0.01 10];
        handles.ax1.YLim=[0.0001 1];
        handles.ax1.XScale='log';
        handles.ax1.YScale='log';
        
    case 'BCHydro2012_5.png'
        handles.e1.String=5.5;
        handles.e2.String=[];
        handles.e3.String=50;
        handles.e4.String=760;
        handles.e5.Value=2;
        handles.e6.Value=1;
        handles.e7.Value=2;
        
        handles.e2.String=50;  plotgmpe(handles);
        handles.e2.String=75;  plotgmpe(handles);
        handles.e2.String=100; plotgmpe(handles);
        handles.e2.String=150; plotgmpe(handles);
        handles.ax1.XLim=[0.01 10];
        handles.ax1.YLim=[0.00001 1];
        handles.ax1.XScale='log';
        handles.ax1.YScale='log';
        
    case 'BCHydro2012_6.png'
        % Figure 11b from Abrahamson 2016
        handles.e1.String=6.5;
        handles.e3.String=50;
        handles.e4.String=760;
        handles.e5.Value=2;
        handles.e6.Value=1;
        handles.e7.Value=2;
        
        handles.e2.String=50;  plotgmpe(handles);
        handles.e2.String=75;  plotgmpe(handles);
        handles.e2.String=100; plotgmpe(handles);
        handles.e2.String=150; plotgmpe(handles);
        handles.ax1.XLim=[0.01 10];
        handles.ax1.YLim=[0.0001 10];
        handles.ax1.XScale='log';
        handles.ax1.YScale='log';
        
    case 'BCHydro2012_7.png'
        % Figure 11c from Abrahamson 2016
        handles.e1.String=7.5;
        handles.e3.String=50;
        handles.e4.String=760;
        handles.e5.Value=2;
        handles.e6.Value=1;
        handles.e7.Value=2;
        
        handles.e2.String=50;  plotgmpe(handles);
        handles.e2.String=75;  plotgmpe(handles);
        handles.e2.String=100; plotgmpe(handles);
        handles.e2.String=150; plotgmpe(handles);
        handles.ax1.XLim=[0.01 10];
        handles.ax1.YLim=[0.0001 10];
        handles.ax1.XScale='log';
        handles.ax1.YScale='log';
end