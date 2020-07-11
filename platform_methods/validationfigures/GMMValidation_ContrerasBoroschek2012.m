function [handles]=GMMValidation_ContrerasBoroschek2012(handles,filename)

handles.t1.String='M';
handles.t2.String='Rrup (km)';
handles.t3.String='Zhyp (km)';
handles.t4.String='Media';

switch filename
    case 'ContrerasBoroschek2012_1.png'
        % Figure 9a Contreras and Boroschek 2012
        handles.e1.String=7.7;
        handles.e2.String=175.6;
        handles.e3.String=38.9;
        handles.e4.Value=2; plotgmpe(handles);
        handles.ax1.XLim=[0 2];
        handles.ax1.YLim=[0 0.429];
        handles.ax1.XScale='linear';
        handles.ax1.YScale='linear';
        
    case 'ContrerasBoroschek2012_2.png'
        % Figure 9b Contreras and Boroschek 2012
        handles.e1.String=7.9;
        handles.e2.String=120.3;
        handles.e3.String=33;
        handles.e4.Value=2; plotgmpe(handles);
        handles.ax1.XLim=[0 2];
        handles.ax1.YLim=[0 0.693];
        handles.ax1.XScale='linear';
        handles.ax1.YScale='linear';
        
    case 'ContrerasBoroschek2012_3.png'
        % Figure 9c Contreras and Boroschek 2012
        handles.e1.String=7.9;
        handles.e2.String=38;
        handles.e3.String=33;
        handles.e4.Value=1; plotgmpe(handles);
        handles.ax1.XLim=[0 2];
        handles.ax1.YLim=[0 0.65];
        handles.ax1.XScale='linear';
        handles.ax1.YScale='linear';
        
    case 'ContrerasBoroschek2012_4.png'
        % Figure 9d Contreras and Boroschek 2012
        handles.e1.String=7.7;
        handles.e2.String=118.5;
        handles.e3.String=38.9;
        handles.e4.Value=1; plotgmpe(handles);
        handles.ax1.XLim=[0 2];
        handles.ax1.YLim=[0 0.3];
        handles.ax1.XScale='linear';
        handles.ax1.YScale='linear';
        
    case 'ContrerasBoroschek2012_5.png'
        % Figure 9d Contreras and Boroschek 2012
        handles.e1.String=8.8;
        handles.e2.String=34;
        handles.e3.String=30;
        handles.e4.Value=2; plotgmpe(handles);
        handles.ax1.XLim=[0 2];
        handles.ax1.YLim=[0 2.29];
        handles.ax1.XScale='linear';
        handles.ax1.YScale='linear';
        
    case 'ContrerasBoroschek2012_6.png'
        % Figure 9d Contreras and Boroschek 2012
        handles.e1.String=8.8;
        handles.e2.String=99;
        handles.e3.String=30;
        handles.e4.Value=2; plotgmpe(handles);
        handles.ax1.XLim=[0 2];
        handles.ax1.YLim=[0 1.578];
        handles.ax1.XScale='linear';
        handles.ax1.YScale='linear';
end