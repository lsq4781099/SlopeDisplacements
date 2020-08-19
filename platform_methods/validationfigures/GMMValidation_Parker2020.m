function [handles]=GMMValidation_Parker2020(handles,filename)

switch filename
    case 'Parker2020_1.png'
        T  = [0.01,0.02,0.03,0.05,0.075,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.75,1,1.5,2,3,4,5,7.5,10];
        lny = nan(4,length(T));
        for i=1:length(T)
            lny(1,i)=Parker2020(T(i),8,30,NaN, 200,0.607,'interface','global');
            lny(2,i)=Parker2020(T(i),8,30,NaN, 400,0.607,'interface','global');
            lny(3,i)=Parker2020(T(i),8,30,NaN, 700,0.607,'interface','global');
            lny(4,i)=Parker2020(T(i),8,30,NaN,1000,0.607,'interface','global');
        end
        plot(handles.ax1,T,exp(lny),'linewidth',1);
        
        handles.ax1.XLim   = [0.01 10];
        handles.ax1.YLim   = [0.001 1];
        handles.ax1.XScale = 'log';
        handles.ax1.YScale = 'log';
        
    case 'Parker2020_2.png'
        T  = [0.01,0.02,0.03,0.05,0.075,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.75,1,1.5,2,3,4,5,7.5,10];
        lny = nan(4,length(T));
        for i=1:length(T)
            lny(1,i)=Parker2020(T(i),7,200,50, 200,0.607,'interface','global');
            lny(2,i)=Parker2020(T(i),7,200,50, 400,0.607,'interface','global');
            lny(3,i)=Parker2020(T(i),7,200,50, 700,0.607,'interface','global');
            lny(4,i)=Parker2020(T(i),7,200,50,1000,0.607,'interface','global');
        end
        plot(handles.ax1,T,exp(lny),'linewidth',1);
        
        handles.ax1.XLim   = [0.01   10];
        handles.ax1.YLim   = [0.0001 0.1];
        handles.ax1.XScale = 'log';
        handles.ax1.YScale = 'log';
        
    case 'Parker2020_3.png'
        nR   = 30;
        Rrup = logsp(25,1000,nR)';
        M    = 9*ones(nR,1);
        lny  = nan(length(Rrup),8);
        
        lny(:,1)=Parker2020(1,M,Rrup,NaN, 760,0.607,'interface','global');
        lny(:,2)=Parker2020(1,M,Rrup,NaN, 760,0.607,'interface','alaska');
        lny(:,3)=Parker2020(1,M,Rrup,NaN, 760,0.607,'interface','aleutian');
        lny(:,4)=Parker2020(1,M,Rrup,NaN, 760,0.607,'interface','central_america_n');
        lny(:,5)=Parker2020(1,M,Rrup,NaN, 760,0.607,'interface','japan_pac');
        lny(:,6)=Parker2020(1,M,Rrup,NaN, 760,0.607,'interface','japan_phi');
        lny(:,7)=Parker2020(1,M,Rrup,NaN, 760,0.607,'interface','south_america_s');
        lny(:,8)=Parker2020(1,M,Rrup,NaN, 760,0.607,'interface','taiwan');
        plot(handles.ax1,Rrup,exp(lny),'linewidth',1);
        
        handles.ax1.XLim   = [25 1000];
        handles.ax1.YLim   = [0.0001 1];
        handles.ax1.XScale = 'log';
        handles.ax1.YScale = 'log';    
        
    case 'Parker2020_4.png'
        T  = [0.01,0.02,0.03,0.05,0.075,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.75,1,1.5,2,3,4,5,7.5,10];
        lny = nan(4,length(T));
        for i=1:length(T)
            lny(1,i)=Parker2020(T(i),8,30,NaN, 200,0.607,'interface','global');
            lny(2,i)=Parker2020(T(i),8,30,NaN, 400,0.607,'interface','global');
            lny(3,i)=Parker2020(T(i),8,30,NaN, 700,0.607,'interface','global');
            lny(4,i)=Parker2020(T(i),8,30,NaN,1000,0.607,'interface','global');
        end
        plot(handles.ax1,T,exp(lny),'linewidth',1);
        
        handles.ax1.XLim   = [0.01 10];
        handles.ax1.YLim   = [0.001 1];
        handles.ax1.XScale = 'log';
        handles.ax1.YScale = 'log';        
        
    case 'Parker2020_5.png'
        T  = [0.01,0.02,0.03,0.05,0.075,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.75,1,1.5,2,3,4,5,7.5,10];
        lny = nan(4,length(T));
        for i=1:length(T)
            lny(1,i)=Parker2020(T(i),7,200,50, 200,0.607,'interface','global');
            lny(2,i)=Parker2020(T(i),7,200,50, 400,0.607,'interface','global');
            lny(3,i)=Parker2020(T(i),7,200,50, 700,0.607,'interface','global');
            lny(4,i)=Parker2020(T(i),7,200,50,1000,0.607,'interface','global');
        end
        plot(handles.ax1,T,exp(lny),'linewidth',1);
        
        handles.ax1.XLim   = [0.01 10];
        handles.ax1.YLim   = [0.0001 1];
        handles.ax1.XScale = 'log';
        handles.ax1.YScale = 'log';        
end