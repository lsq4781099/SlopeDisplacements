function [handles]=GMMValidation_Kanno2006(handles,filename)


switch filename
    case 'Kanno2006_1.png'
        T    = 0; %PGA
        Rrup = logsp(1,300,100);
        on   = ones(1,100);
        Zhyp = 5*on;
        lny  = zeros(4,100);
        lny(1,:) = Kanno2006(T,6.5,Rrup,Zhyp,300);
        lny(2,:) = Kanno2006(T,7.0,Rrup,Zhyp,300);
        lny(3,:) = Kanno2006(T,7.5,Rrup,Zhyp,300);
        lny(4,:) = Kanno2006(T,8.0,Rrup,Zhyp,300);
        
        PGA = exp(lny);
        plot(handles.ax1,Rrup,PGA*980.66,'b-','linewidth',1,'tag','curves')
        set(handles.ax1,'xscale','log','yscale','log','xlim',[1 300],'ylim',[1e-1 1e4])
        xlabel(handles.ax1,'Fault Distance (km)')
        ylabel(handles.ax1,'PGA (cm/s2)')
        
    case 'Kanno2006_2.png'
        T    = -1; %PGV
        Rrup = logsp(1,300,100);
        on   = ones(1,100);
        Zhyp = 5*on;
        lny  = zeros(4,100);
        lny(1,:) = Kanno2006(T,6.5,Rrup,Zhyp,300);
        lny(2,:) = Kanno2006(T,7.0,Rrup,Zhyp,300);
        lny(3,:) = Kanno2006(T,7.5,Rrup,Zhyp,300);
        lny(4,:) = Kanno2006(T,8.0,Rrup,Zhyp,300);
        
        PGA = exp(lny);
        plot(handles.ax1,Rrup,PGA,'b-','linewidth',1,'tag','curves')
        set(handles.ax1,'xscale','log','yscale','log','xlim',[1 300],'ylim',[1e-1 1e4])
        xlabel(handles.ax1,'Fault Distance (km)')
        ylabel(handles.ax1,'PVA (cm/s)')
        
        
    case 'Kanno2006_3.png'
        T    = 0; %PGA
        Rrup = logsp(30,300,100);
        on   = ones(1,100);
        Zhyp = 50*on;
        lny  = zeros(4,100);
        lny(1,:) = Kanno2006(T,6.5,Rrup,Zhyp,300);
        lny(2,:) = Kanno2006(T,7.0,Rrup,Zhyp,300);
        lny(3,:) = Kanno2006(T,7.5,Rrup,Zhyp,300);
        lny(4,:) = Kanno2006(T,8.0,Rrup,Zhyp,300);
        
        PGA = exp(lny);
        plot(handles.ax1,Rrup,PGA*980.66,'b-','linewidth',1,'tag','curves')
        set(handles.ax1,'xscale','log','yscale','log','xlim',[10 300],'ylim',[1e-1 1e4])
        xlabel(handles.ax1,'Fault Distance (km)')
        ylabel(handles.ax1,'PGA (cm/s2)')
        
    case 'Kanno2006_4.png'
        T    = -1; %PGV
        Rrup = logsp(30,300,100);
        on   = ones(1,100);
        Zhyp = 50*on;
        lny  = zeros(4,100);
        lny(1,:) = Kanno2006(T,6.5,Rrup,Zhyp,300);
        lny(2,:) = Kanno2006(T,7.0,Rrup,Zhyp,300);
        lny(3,:) = Kanno2006(T,7.5,Rrup,Zhyp,300);
        lny(4,:) = Kanno2006(T,8.0,Rrup,Zhyp,300);
        
        PGA = exp(lny);
        plot(handles.ax1,Rrup,PGA,'b-','linewidth',1,'tag','curves')
        set(handles.ax1,'xscale','log','yscale','log','xlim',[10 300],'ylim',[1e-1 1e4])
        xlabel(handles.ax1,'Fault Distance (km)')
        ylabel(handles.ax1,'PVA (cm/s)')
        
    case 'Kanno2006_5.png'
        handles.e1.String=7;
        handles.e4.String=300;
        handles.e2.String=5;  handles.e3.String=5; plotgmpe(handles,980.66);
        handles.e2.String=20; handles.e3.String=5; plotgmpe(handles,980.66);
        handles.e2.String=40; handles.e3.String=5; plotgmpe(handles,980.66);
        handles.e2.String=80; handles.e3.String=5; plotgmpe(handles,980.66);
        
        handles.ax1.XLim=[0.05  5];
        handles.ax1.YLim=[1e-1 1e4];
        handles.ax1.XScale='log';
        handles.ax1.YScale='log';
        
    case 'Kanno2006_6.png'
        handles.e1.String=7;
        handles.e4.String=300;
        handles.e2.String=40;  handles.e3.String=50; plotgmpe(handles,980.66);
        handles.e2.String=80;  handles.e3.String=50; plotgmpe(handles,980.66);
        handles.e2.String=120; handles.e3.String=50; plotgmpe(handles,980.66);
        handles.e2.String=200; handles.e3.String=50; plotgmpe(handles,980.66);
        
        handles.ax1.XLim=[0.05  5];
        handles.ax1.YLim=[1e-1 1e4];
        handles.ax1.XScale='log';
        handles.ax1.YScale='log';
        
end