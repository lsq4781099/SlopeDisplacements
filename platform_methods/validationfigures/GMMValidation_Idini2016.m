function [handles]=GMMValidation_Idini2016(handles,filename)

switch filename
    case 'Idini2016_1.png'
        % Figure 9a from Idini 2016
        handles.e1.String=8.5;
        handles.e2.String=50;
        handles.e3.String=[];
        handles.e4.String=400;
        handles.e5.Value=1;
        handles.e6.Value=1; plotgmpe(handles);
        handles.e6.Value=2; plotgmpe(handles);
        handles.e6.Value=5; plotgmpe(handles);
        handles.ax1.XLim=[0.01 10];
        handles.ax1.YLim=[0.001 10];
        handles.ax1.XScale='log';
        handles.ax1.YScale='log';
        
    case 'Idini2016_2.png'
        % Figure 9b from Idini 2016
        handles.e1.String=7.5;
        handles.e2.String=100;
        handles.e3.String=100;
        handles.e4.String=400;
        handles.e5.Value=2;
        handles.e6.Value=1; plotgmpe(handles);
        handles.e6.Value=2; plotgmpe(handles);
        handles.e6.Value=5; plotgmpe(handles);        

        handles.ax1.XLim=[0.01 10];
        handles.ax1.YLim=[0.001 10];
        handles.ax1.XScale='log';
        handles.ax1.YScale='log';
        
    case 'Idini2016_3.png'
        % Figure 10a from Idini 2016
        handles.e2.String=30;
        handles.e3.String=[];
        handles.e4.String=1100;
        handles.e5.Value=1;
        handles.e6.Value=1;
        
        handles.e1.String=9; plotgmpe(handles);
        handles.e1.String=8; plotgmpe(handles);
        handles.e1.String=7; plotgmpe(handles);
        handles.e1.String=6; plotgmpe(handles);
        handles.ax1.XLim=[0.01 10];
        handles.ax1.YLim=[0.001 5];
        handles.ax1.XScale='log';
        handles.ax1.YScale='log';
        
    case 'Idini2016_4.png'
        % Figure 10b from Idini 2016
        handles.e2.String=60;
        handles.e3.String=60;
        handles.e4.String=1100;
        handles.e5.Value=2;
        handles.e6.Value=1;
        
        handles.e1.String=8; plotgmpe(handles);
        handles.e1.String=7; plotgmpe(handles);
        handles.e1.String=6; plotgmpe(handles);
    
    case 'Idini2016_5.png'
        % Figure 10c from Idini 2016
        handles.e2.String=100;
        handles.e3.String=[];
        handles.e4.String=1100;
        handles.e5.Value=1;
        handles.e6.Value=1;
        
        handles.e1.String=9; plotgmpe(handles);
        handles.e1.String=8; plotgmpe(handles);
        handles.e1.String=7; plotgmpe(handles);
        handles.e1.String=6; plotgmpe(handles);
        handles.ax1.XLim=[0.01 10];
        handles.ax1.YLim=[0.001 5];
        handles.ax1.XScale='log';
        handles.ax1.YScale='log';
        
        
    case 'Idini2016_6.png'
        % Figure 10b from Idini 2016
        handles.e2.String=100;
        handles.e3.String=100;
        handles.e4.String=1100;
        handles.e5.Value=2;
        handles.e6.Value=1;
        
        handles.e1.String=8; plotgmpe(handles);
        handles.e1.String=7; plotgmpe(handles);
        handles.e1.String=6; plotgmpe(handles);
%         
end