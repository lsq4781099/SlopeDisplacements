function [handles]=GMMValidation_I2008(handles,filename)

switch filename

    case 'I2008_1.png'
        handles.e1.String=7;
        handles.e2.String=6;
        handles.e3.String='500';
        handles.e4.Value=4;
        plotgmpe(handles);
        handles.ax1.XLim=[0.01  10];
        handles.ax1.YLim=[0.001 10];
        handles.ax1.XScale='log';
        handles.ax1.YScale='log';    
    
    case 'I2008_2.png'
        handles.e1.String=6;
        handles.e2.String=6;
        handles.e3.String='500'; 
        handles.e4.Value=1;
        
        plotgmpe(handles);
        handles.ax1.XLim=[0.01  10];
        handles.ax1.YLim=[0.001 10];
        handles.ax1.XScale='log';
        handles.ax1.YScale='log';
end