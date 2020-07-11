function [handles]=GMMValidation_AS2008(handles,filename)

switch filename
    case 'AS2008_1.png'
        %1  2     3    4   5    6     7   8   9    10         11     12
        %M, rrup, rjb, rx, dip, ztor, W, Z10, Vs30,mechanism, event, Vs30type
        handles.e3.String = 30; %rjb
        handles.e4.String = 30; %rx
        handles.e6.String = 90; %dip
        handles.e7.String = 10;
        handles.e8.String = 270;
        handles.e10.Value = 1; %{'strike-slip','normal','normal-oblique','reverse','reverse-oblique','thrust'};
        handles.e11.Value = 2; %{'aftershock','mainshock','foreshock','swarms'};
        handles.e12.Value = 1; % {'measured','inferred'};
        handles.e1.String = 5; handles.e2.String = 30; handles.e5.String = 7; handles.e9.String = Z10_default_AS08_NGA(270); plotgmpe(handles);
        handles.e1.String = 6; handles.e2.String = 30; handles.e5.String = 3; handles.e9.String = Z10_default_AS08_NGA(270); plotgmpe(handles);
        handles.e1.String = 7; handles.e2.String = 30; handles.e5.String = 0; handles.e9.String = Z10_default_AS08_NGA(270); plotgmpe(handles);
        handles.e1.String = 8; handles.e2.String = 30; handles.e5.String = 0; handles.e9.String = Z10_default_AS08_NGA(270); plotgmpe(handles);
        
        handles.ax1.XLim=[0.01  10];
        handles.ax1.YLim=[0.0001 1];
        handles.ax1.XScale='log';
        handles.ax1.YScale='log';
        
    case 'AS2008_2.png'
        %1-M, 2-Vs30, 3-Rrup, 4-Rjb, 5-Rx, 6-dip, 7-Ztor, 8-Z10, 9-W, 10-mechanism, 11-eventtype, 12-side, 13-Vs30source
        handles.e3.String = 30;
        handles.e4.String = 30;
        handles.e6.String = 90;
        handles.e7.String = 10;
        handles.e8.String = 550;
        handles.e10.Value = 1; %{'strike-slip','normal','normal-oblique','reverse','reverse-oblique','thrust'};
        handles.e11.Value = 2; %{'aftershock','mainshock','foreshock','swarms'};
        handles.e12.Value = 1; %{'hangingwall','foot wall','other'};
        handles.e13.Value = 1; % {'measured','inferred'};
        handles.e1.String=5; handles.e2.String = 30; handles.e5.String = 7; handles.e9.String = Z10_default_AS08_NGA(550); plotgmpe(handles);
        handles.e1.String=6; handles.e2.String = 30; handles.e5.String = 3; handles.e9.String = Z10_default_AS08_NGA(550); plotgmpe(handles);
        handles.e1.String=7; handles.e2.String = 30; handles.e5.String = 0; handles.e9.String = Z10_default_AS08_NGA(550); plotgmpe(handles);
        handles.e1.String=8; handles.e2.String = 30; handles.e5.String = 0; handles.e9.String = Z10_default_AS08_NGA(550); plotgmpe(handles);
        
        handles.ax1.XLim=[0.01  10];
        handles.ax1.YLim=[1e-5 1];
        handles.ax1.XScale='log';
        handles.ax1.YScale='log';
        
    case 'AS2008_3.png'
        handles.e3.String = 30;
        handles.e4.String = 30;
        handles.e6.String = 90;
        handles.e7.String = 10;
        handles.e8.String = 760;
        handles.e10.Value = 1; %{'strike-slip','normal','normal-oblique','reverse','reverse-oblique','thrust'};
        handles.e11.Value = 2; %{'aftershock','mainshock','foreshock','swarms'};
        handles.e12.Value = 1; %{'hangingwall','foot wall','other'};
        handles.e13.Value = 1; % {'measured','inferred'};
        handles.e1.String=5; handles.e2.String = 30; handles.e5.String = 7; handles.e9.String = Z10_default_AS08_NGA(760); plotgmpe(handles);
        handles.e1.String=6; handles.e2.String = 30; handles.e5.String = 3; handles.e9.String = Z10_default_AS08_NGA(760); plotgmpe(handles);
        handles.e1.String=7; handles.e2.String = 30; handles.e5.String = 0; handles.e9.String = Z10_default_AS08_NGA(760); plotgmpe(handles);
        handles.e1.String=8; handles.e2.String = 30; handles.e5.String = 0; handles.e9.String = Z10_default_AS08_NGA(760); plotgmpe(handles);
        
        handles.ax1.XLim=[0.01  10];
        handles.ax1.YLim=[1e-5 1];
        handles.ax1.XScale='log';
        handles.ax1.YScale='log';
        
    case 'AS2008_4.png'
        handles.e3.String = 1;
        handles.e4.String = 1;
        handles.e6.String = 90;
        handles.e7.String = 15;
        handles.e8.String = 270;
        handles.e10.Value = 1; %{'strike-slip','normal','normal-oblique','reverse','reverse-oblique','thrust'};
        handles.e11.Value = 2; %{'aftershock','mainshock','foreshock','swarms'};
        handles.e12.Value = 1; %{'hangingwall','foot wall','other'};
        handles.e13.Value = 1; % {'measured','inferred'};
        handles.e1.String=5; handles.e2.String = 1; handles.e5.String = 7; handles.e9.String = Z10_default_AS08_NGA(270); plotgmpe(handles);
        handles.e1.String=6; handles.e2.String = 1; handles.e5.String = 3; handles.e9.String = Z10_default_AS08_NGA(270); plotgmpe(handles);
        handles.e1.String=7; handles.e2.String = 1; handles.e5.String = 0; handles.e9.String = Z10_default_AS08_NGA(270); plotgmpe(handles);
        handles.e1.String=8; handles.e2.String = 1; handles.e5.String = 0; handles.e9.String = Z10_default_AS08_NGA(270); plotgmpe(handles);
        
        handles.ax1.XLim=[0.01  10];
        handles.ax1.YLim=[1e-4 10];
        handles.ax1.XScale='log';
        handles.ax1.YScale='log';
        
    case 'AS2008_5.png'
        %1-M, 2-Vs30, 3-Rrup, 4-Rjb, 5-Rx, 6-dip, 7-Ztor, 8-Z10, 9-W, 10-mechanism, 11-eventtype, 12-side, 13-Vs30source
        
        handles.e3.String = 1;
        handles.e4.String = 1;
        handles.e6.String = 90;
        handles.e7.String = 15;
        handles.e8.String = 550;
        handles.e10.Value = 1; %{'strike-slip','normal','normal-oblique','reverse','reverse-oblique','thrust'};
        handles.e11.Value = 2; %{'aftershock','mainshock','foreshock','swarms'};
        handles.e12.Value = 1; %{'hangingwall','foot wall','other'};
        handles.e13.Value = 1; % {'measured','inferred'};
        handles.e1.String=5; handles.e2.String = 1; handles.e5.String = 7; handles.e9.String = Z10_default_AS08_NGA(550); plotgmpe(handles);
        handles.e1.String=6; handles.e2.String = 1; handles.e5.String = 3; handles.e9.String = Z10_default_AS08_NGA(550); plotgmpe(handles);
        handles.e1.String=7; handles.e2.String = 1; handles.e5.String = 0; handles.e9.String = Z10_default_AS08_NGA(550); plotgmpe(handles);
        handles.e1.String=8; handles.e2.String = 1; handles.e5.String = 0; handles.e9.String = Z10_default_AS08_NGA(550); plotgmpe(handles);
        
        handles.ax1.XLim=[0.01  10];
        handles.ax1.YLim=[1e-4 10];
        handles.ax1.XScale='log';
        handles.ax1.YScale='log';
        
end