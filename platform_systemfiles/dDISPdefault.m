function handles=dDISPdefault(handles,txt,edi)

handles.e1.BackgroundColor=[1 1 1];
handles.e2.BackgroundColor=[1 1 1];
handles.e3.BackgroundColor=[1 1 1];
handles.e4.BackgroundColor=[1 1 1];
handles.e5.BackgroundColor=[1 1 1];
handles.e6.BackgroundColor=[1 1 1];
handles.e7.BackgroundColor=[1 1 1];
handles.e8.BackgroundColor=[1 1 1];
handles.e9.BackgroundColor=[1 1 1];
handles.e10.BackgroundColor=[1 1 1];
handles.e11.BackgroundColor=[1 1 1];

%% SUBDUCTION EARTHQUAKES
meth = handles.meth;
val  = handles.Dpop.Value;
str  = meth(val).label;
handles.fun = meth(val).func;
switch str
    case 'BMT 2017 Sa(M)'
        set(txt(1:2),'Visible','on');
        set(edi(1:2),'Visible','on');
        handles.t1.String='Sa(1.5Ts)';
        handles.t2.String='Magnitude';
        handles.e1.String='0.5';
        handles.e2.String='7.0';
        handles.e1.BackgroundColor=[1 1 0.7];
        handles.e2.BackgroundColor=[1 1 0.7];
        xlabel(handles.ax1,'Sa(1.5Ts) (g)')
        handles.im_1 = 1e-2;
        handles.im_2 = 1.0;
    case 'BT 2007 Sa'
        set(txt(1:1),'Visible','on');
        set(edi(1:1),'Visible','on');
        handles.t1.String='Sa(1.5Ts)';
        handles.e1.String='0.5';
        handles.e1.BackgroundColor=[1 1 0.7];
        xlabel(handles.ax1,'Sa(1.5Ts) (g)')
        handles.im_1 = 1e-2;
        handles.im_2 = 1.0;
    case 'BT 2007 Sa(M)'
        set(txt(1:2),'Visible','on');
        set(edi(1:2),'Visible','on');
        handles.t1.String='Sa(1.5Ts)';
        handles.t2.String='Magnitude';
        handles.e1.String='0.5';
        handles.e2.String='7.0';
        handles.e1.BackgroundColor=[1 1 0.7];
        handles.e2.BackgroundColor=[1 1 0.7];
        xlabel(handles.ax1,'Sa(1.5Ts) (g)')
        handles.im_1 = 1e-2;
        handles.im_2 = 1.0;
    case 'BM 2019 NonNF (M)'
        set(txt(1:2),'Visible','on');
        set(edi(1:2),'Visible','on');
        handles.t1.String='Sa(1.3Ts)';
        handles.t2.String='Magnitude';
        handles.e1.String='0.5';
        handles.e2.String='7.0';
        handles.e1.BackgroundColor=[1 1 0.7];
        handles.e2.BackgroundColor=[1 1 0.7];
        xlabel(handles.ax1,'Sa(1.3Ts) (g)')
        handles.im_1 = 1e-2;
        handles.im_2 = 1.0;
    case 'Jibson  2007 (M)'
        set(txt(1:2),'Visible','on');
        set(edi(1:2),'Visible','on');
        handles.t1.String='PGA (g)';
        handles.t2.String='Magnitude';
        handles.e1.String='0.5';
        handles.e2.String='7.0';
        handles.e1.BackgroundColor=[1 1 0.7];
        handles.e2.BackgroundColor=[1 1 0.7];
        xlabel(handles.ax1,'PGA (g)')
        handles.im_1 = 1e-2;
        handles.im_2 = 1.0;
    case 'Jibson  2007 Ia'
        set(txt(1:1),'Visible','on');
        set(edi(1:1),'Visible','on');
        handles.t1.String='Ia (cm/s)';
        handles.e1.String='500';
        handles.e1.BackgroundColor=[1 1 0.7];
        xlabel(handles.ax1,'Ia (cm/s)')
        handles.im_1 = 1;
        handles.im_2 = 1000;
    case 'RA 2011 (Rigid)'
        set(txt(1:2),'Visible','on');
        set(edi(1:2),'Visible','on');
        handles.t1.String='PGV (cm/s)';
        handles.t2.String='PGA (g)';
        handles.e1.String='60';
        handles.e2.String='0.5';
        handles.e1.BackgroundColor=[1 1 0.7];
        handles.e2.BackgroundColor=[1 1 0.7];
        xlabel(handles.ax1,'PGA (g)')
        handles.im_1 = 1e-2;
        handles.im_2 = 1.0;
    case 'RA 2011 (Flexible)'
        set(txt(1:2),'Visible','on');
        set(edi(1:2),'Visible','on');
        handles.t1.String='kvmax (cm/s)';
        handles.t2.String='kmax (g)';
        handles.e1.String='60';
        handles.e2.String='0.5';
        handles.e1.BackgroundColor=[1 1 0.7];
        handles.e2.BackgroundColor=[1 1 0.7];
        xlabel(handles.ax1,'kmax (g)')
        handles.im_1 = 1e-2;
        handles.im_2 = 1.0;
    case 'RS 2009 (Scalar-M)'
        set(txt(1:2),'Visible','on');
        set(edi(1:2),'Visible','on');
        handles.t1.String='PGA (g)';
        handles.t2.String='Magnitude';
        handles.e1.String='0.5';
        handles.e2.String='7.0';
        handles.e1.BackgroundColor=[1 1 0.7];
        handles.e2.BackgroundColor=[1 1 0.7];
        xlabel(handles.ax1,'PGA (g)')
        handles.im_1 = 1e-2;
        handles.im_2 = 1.0;
    case 'RS 2009 (Vector)'
        set(txt(1:2),'Visible','on');
        set(edi(1:2),'Visible','on');
        handles.t1.String='PGV (cm/s)';
        handles.t2.String='PGA (g)';
        handles.e1.String='60';
        handles.e2.String='0.5';
        handles.e1.BackgroundColor=[1 1 0.7];
        handles.e2.BackgroundColor=[1 1 0.7];
        xlabel(handles.ax1,'PGA (g)')
        handles.im_1 = 1e-2;
        handles.im_2 = 1.0;
    case 'AM 1988'
        set(txt(1:1),'Visible','on');
        set(edi(1:1),'Visible','on');
        handles.t1.String='PGA (g)';
        handles.e1.String='0.5';
        handles.e1.BackgroundColor=[1 1 0.7];
        xlabel(handles.ax1,'PGA (g)')
        handles.im_1 = 1e-2;
        handles.im_2 = 1.0;

end

