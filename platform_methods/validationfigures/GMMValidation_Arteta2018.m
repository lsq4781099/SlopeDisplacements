function [handles]=GMMValidation_Arteta2018(handles,filename)

switch filename
    case 'Arteta2018_1.png'
        % Figure 10a from Abrahamson 2016
        handles.e1.String=7.2;
        handles.e2.String=120;
        handles.e3.Value=1;
        handles.epsilon = [-1 0 1];
        plotgmpe(handles);
        handles.ax1.XLim=[0.01 10];
        handles.ax1.YLim=[0.0001 1];
        handles.ax1.XScale='log';
        handles.ax1.YScale='log';
        xlabel(handles.ax1,'Period [s]')
        ylabel(handles.ax1,'Sa RotD50s [g]')
        
    case  'Arteta2018_2.png'
        handles.e2.String=80;
        handles.e3.Value=1;
        handles.e1.String=5; plotgmpe(handles);
        handles.e1.String=6; plotgmpe(handles);
        handles.e1.String=7; plotgmpe(handles);
        handles.e1.String=8; plotgmpe(handles);
        
        handles.ax1.ColorOrderIndex=1;
        handles.e3.Value=2;
        handles.e1.String=5; plotgmpe(handles);
        handles.e1.String=6; plotgmpe(handles);
        handles.e1.String=7; plotgmpe(handles);
        handles.e1.String=8; plotgmpe(handles);
        
        handles.ax1.XLim=[0.01 10];
        handles.ax1.YLim=[0.00001 1];
        handles.ax1.XScale='log';
        handles.ax1.YScale='log';
        xlabel(handles.ax1,'Period [s]')
        ylabel(handles.ax1,'Sa RotD50s [g]')
        
    case  'Arteta2018_3.png'
        handles.e2.String=120;
        handles.e3.Value=1;
        handles.e1.String=5; plotgmpe(handles);
        handles.e1.String=6; plotgmpe(handles);
        handles.e1.String=7; plotgmpe(handles);
        handles.e1.String=8; plotgmpe(handles);
        
        handles.ax1.ColorOrderIndex=1;
        handles.e3.Value=2;
        handles.e1.String=5; plotgmpe(handles);
        handles.e1.String=6; plotgmpe(handles);
        handles.e1.String=7; plotgmpe(handles);
        handles.e1.String=8; plotgmpe(handles);
        
        handles.ax1.XLim=[0.01 10];
        handles.ax1.YLim=[0.00001 1];
        handles.ax1.XScale='log';
        handles.ax1.YScale='log';
        xlabel(handles.ax1,'Period [s]')
        ylabel(handles.ax1,'Sa RotD50s [g]')
        
    case  'Arteta2018_4.png'
        handles.e2.String=200;
        handles.e3.Value=1;
        handles.e1.String=5; plotgmpe(handles);
        handles.e1.String=6; plotgmpe(handles);
        handles.e1.String=7; plotgmpe(handles);
        handles.e1.String=8; plotgmpe(handles);
        
        handles.ax1.ColorOrderIndex=1;
        handles.e3.Value=2;
        handles.e1.String=5; plotgmpe(handles);
        handles.e1.String=6; plotgmpe(handles);
        handles.e1.String=7; plotgmpe(handles);
        handles.e1.String=8; plotgmpe(handles);
        
        handles.ax1.XLim=[0.01 10];
        handles.ax1.YLim=[0.00001 1];
        handles.ax1.XScale='log';
        handles.ax1.YScale='log'; 
        xlabel(handles.ax1,'Period [s]')
        ylabel(handles.ax1,'Sa RotD50s [g]')        
        
        
    case  'Arteta2018_5.png'
        M    = linspace(5,8,50)';
        on   = ones(50,1);
        lny(:,1)=Arteta2018(0,M,80*on,'rock','forearc');
        lny(:,2)=Arteta2018(0,M,160*on,'rock','forearc');
        lny(:,3)=Arteta2018(0,M,300*on,'rock','forearc');
        lny(:,4)=Arteta2018(0,M,500*on,'rock','forearc');
        
        plot(handles.ax1,M,exp(lny),'linewidth',1)
        handles.ax1.XLim=[5 8];
        handles.ax1.YLim=[0.0001 0.1];
        handles.ax1.XScale='linear';
        handles.ax1.YScale='log';
        xlabel(handles.ax1,'Magnitude, Mw')
        ylabel(handles.ax1,'PGA [g]')
        
    case 'Arteta2018_6.png'        
        handles.e1.String=6;
        handles.e3.Value=1;
        
        handles.e4.Value=1;
        handles.e2.String=80;  plotgmpe(handles);
        handles.e2.String=160; plotgmpe(handles);
        handles.e2.String=300; plotgmpe(handles);
        handles.e2.String=500; plotgmpe(handles);
        
        handles.ax1.ColorOrderIndex=1;
        handles.e4.Value=2;
        handles.e2.String=80;  plotgmpe(handles);
        handles.e2.String=160; plotgmpe(handles);
        handles.e2.String=300; plotgmpe(handles);
        handles.e2.String=500; plotgmpe(handles);        
        
        handles.ax1.XLim=[0.01 10];
        handles.ax1.YLim=[0.00001 1];
        handles.ax1.XScale='log';
        handles.ax1.YScale='log'; 
        xlabel(handles.ax1,'Period [s]')
        ylabel(handles.ax1,'Sa RotD50s [g]')    
        
    case  'Arteta2018_7.png'
        Rhyp  = logsp(100,500,10)';
        on    = ones(10,1);
        [lny(:,1),sig(:,1)]=Arteta2018(0,6.3*on,Rhyp,'rock','forearc');
        [lny(:,2),sig(:,2)]=Arteta2018(0,6.3*on,Rhyp,'soil','forearc');
        
        plot(handles.ax1,Rhyp,exp(lny),'linewidth',1)
        handles.ax1.ColorOrderIndex=1;
        plot(handles.ax1,Rhyp,exp(lny+1.5*sig),'--','linewidth',1)
        handles.ax1.ColorOrderIndex=1;
        plot(handles.ax1,Rhyp,exp(lny-1.5*sig),'--','linewidth',1)
        handles.ax1.XLim=[100 1000];
        handles.ax1.YLim=[0.001 1];
        handles.ax1.XScale='log';
        handles.ax1.YScale='log';
        xlabel(handles.ax1,'Rhyp [km]')
        ylabel(handles.ax1,'PGA [g]')        
        
        
        
end