function [handles]=GMMValidation_MontalvaBastias2017(handles,filename)

switch filename
    case 'MontalvaBastias2017_1.png'
        % Figure 3a from Montalva and Bastias 2017
        handles.e2.String=25;
        handles.e3.String=25;
        handles.e4.String=300;
        handles.e5.Value=1;
        handles.e6.Value=1;
        
        handles.e1.String=6.5; plotgmpe(handles);
        handles.e1.String=8.5; plotgmpe(handles);
        handles.ax1.XLim=[0.01 10];
        handles.ax1.YLim=[0.001 10];
        handles.ax1.XScale='log';
        handles.ax1.YScale='log';
        
    case 'MontalvaBastias2017_2.png'
        % Figure 3b from Montalva and Bastias 2017
        handles.e2.String=50;
        handles.e3.String=50;
        handles.e4.String=300;
        handles.e5.Value=1;
        handles.e6.Value=1;
        
        handles.e1.String=6.5; plotgmpe(handles);
        handles.e1.String=8.5; plotgmpe(handles);
        handles.ax1.XLim=[0.01 10];
        handles.ax1.YLim=[0.001 10];
        handles.ax1.XScale='log';
        handles.ax1.YScale='log';
        
    case 'MontalvaBastias2017_3.png'
        % Figure 3c from Montalva and Bastias 2017
        handles.e2.String=100;
        handles.e3.String=100;
        handles.e4.String=300;
        handles.e5.Value=1;
        handles.e6.Value=1;
        
        handles.e1.String=6.5; plotgmpe(handles);
        handles.e1.String=8.5; plotgmpe(handles);
        handles.ax1.XLim=[0.01 10];
        handles.ax1.YLim=[0.0001 1];
        handles.ax1.XScale='log';
        handles.ax1.YScale='log';
        
    case 'MontalvaBastias2017_4.png'
        % Figure 3d from Montalva and Bastias 2017
        handles.e2.String=150;
        handles.e3.String=150;
        handles.e4.String=300;
        handles.e5.Value=1;
        handles.e6.Value=1;
        
        handles.e1.String=6.5; plotgmpe(handles);
        handles.e1.String=8.5; plotgmpe(handles);
        handles.ax1.XLim=[0.01 10];
        handles.ax1.YLim=[0.0001 1];
        handles.ax1.XScale='log';
        handles.ax1.YScale='log';
        
    case 'MontalvaBastias2017_5.png'
        % Figure 4a from Montalva and Bastias 2017
        handles.e2.String=75;
        handles.e3.String=100;
        handles.e4.String=300;
        handles.e5.Value=2;
        handles.e6.Value=1;
        
        handles.e1.String=6.5; plotgmpe(handles);
        handles.e1.String=8.5; plotgmpe(handles);
        handles.ax1.XLim=[0.01 10];
        handles.ax1.YLim=[0.0001 10];
        handles.ax1.XScale='log';
        handles.ax1.YScale='log';
        
    case 'MontalvaBastias2017_6.png'
        % Figure 4b from Montalva and Bastias 2017
        handles.e2.String=100;
        handles.e3.String=100;
        handles.e4.String=300;
        handles.e5.Value=2;
        handles.e6.Value=1;
        
        handles.e1.String=6.5; plotgmpe(handles);
        handles.e1.String=8.5; plotgmpe(handles);
        handles.ax1.XLim=[0.01 10];
        handles.ax1.YLim=[0.0001 10];
        handles.ax1.XScale='log';
        handles.ax1.YScale='log';
        
    case 'MontalvaBastias2017_7.png'
        % Figure 4c from Montalva and Bastias 2017
        handles.e2.String=150;
        handles.e3.String=100;
        handles.e4.String=300;
        handles.e5.Value=2;
        handles.e6.Value=1;
        
        handles.e1.String=6.5; plotgmpe(handles);
        handles.e1.String=8.5; plotgmpe(handles);
        handles.ax1.XLim=[0.01 10];
        handles.ax1.YLim=[1e-5 1];
        handles.ax1.XScale='log';
        handles.ax1.YScale='log';
        
    case 'MontalvaBastias2017_8.png'
        % Figure 4d from Montalva and Bastias 2017
        handles.e2.String=200;
        handles.e3.String=100;
        handles.e4.String=300;
        handles.e5.Value=2;
        handles.e6.Value=1;
        
        handles.e1.String=6.5; plotgmpe(handles);
        handles.e1.String=8.5; plotgmpe(handles);
        handles.ax1.XLim=[0.01 10];
        handles.ax1.YLim=[1e-5 1];
        handles.ax1.XScale='log';
        handles.ax1.YScale='log';
end