function [LBS,HL_BI,HL_R15] = calc_LBS_FS (Data, wt, Df,Mw, PGA)
%% SOIL LIQUEFACTION TRIGGERING EVALUATION
% This MATLAB code estimates liquefaction triggering based on Boulanger &
% Idriss (2014) method. It also estimates Relative Density according to
% Boualanger ad Idriss (2008), Kulhawy and Mayne (1990)
% and Jamiolkowsky (2001)

% Create a Microsoft Excel File (*.xls) with:
%       * Depth in "m" in column "A"
%       * qc in "MPa" in column "B"
%       * fs in "MPa" in column "C"
%       * u2 in "MPa" in column "D"`
% Make the name of the excel file, the name of the CPT

%% Input Parameters
d_in        = Data(:,1);
qc_in       = Data(:,2);
fs_in       = Data(:,3);
u2_in       = Data(:,4);

a           = 0.58;            % Net Area Ratio
uw_awt      = 18;              % Assume Unit Weight above WT
cfc         = 0.07;            % Cfc for site specific Ic-FC relationship (0.07 for chch)
p_l         = 0.15;            % Use 0.15 for desing or 0.5 for back analysis
NA          = numel(wt);       % Number of Analyses
Pa          = 101.3;           % Atmospheric pressure in kPa
g_w         = 9.81;            % Unit Weight of water
d_d         = d_in(2)-d_in(1); % Spacing in depth

%% Remove extra points becuase shaft cell is X cm above tip cell
N_d   = numel(d_in);
N_qc  = numel(qc_in);
N_fs  = numel(fs_in);
N_u2  = numel(u2_in);

N1 = min(N_d,N_qc);
N2 = min(N1,N_fs);
N3  = min(N2,N_u2);

d_2   = d_in(1:N3);
qc_2  = qc_in(1:N3);
fs_2  = fs_in(1:N3);
u2_2  = u2_in(1:N3);

%% Add points above the first measured point if neccessary
if d_2(1)>=d_d
    N_ex = fix((d_2(1)/d_d)-1);
else
    N_ex = 0;
end

N = N_ex + N3;    % Number of rows in each vector for the rest of analysis
d   = (1:N)'*d_d;
qc  = [zeros(N_ex,1);qc_2]; qc(qc==0)=0.001;
fs  = [zeros(N_ex,1);fs_2]; fs(fs==0)=0.0001;
u2  = [zeros(N_ex,1);u2_2]; u2(u2==0)=0.00001;

%% Calculate state of stresses
% Unit Weight from Robertson (2010)
Rf          = nan([N,1]);            % Friction Ratio (%)
qt          = nan([N,1]);            % Corrected tip resistance
g_r         = uw_awt*ones([N,1]);            % Unit Weight

qt(N_ex+1:N)  = 1000*qc(N_ex+1:N) + 1000*u2(N_ex+1:N)*(1-a);
Rf(N_ex+1:N)  = (1000*fs(N_ex+1:N)./qt(N_ex+1:N))*100;
g_r(N_ex+1:N) = g_w*(0.27*log10(Rf(N_ex+1:N))+0.36*log10(qt(N_ex+1:N)/Pa)+1.236);

%Svo, u and S'vo
svo     = zeros([N,1]);
u       = zeros([N,NA]);
s_vo    = zeros([N,NA]);

for j = 1:NA
    for i = 1:N
        if i == 1
            svo(i) = g_r(i)*d(i);
        else
            svo(i) = svo(i-1)+d_d*g_r(i);
        end
        if d(i) <= wt(j,1)
            u(i,j)    = 0;
        else
            u(i,j)    = g_w*(d(i)-wt(j,1));
        end
        s_vo(i,j) = svo(i) - u(i,j);
    end
end

%% Estimation of Ic and Robertson Parameters
c_n     = zeros([N,NA]);
Qtn     = zeros([N,NA]);
Fr      = zeros([N,NA]);
ic      = zeros([N,NA]);
n_exp   = zeros([N,NA]);
delta_n = zeros([N,NA])+1;
toler   = 0.001;
n_ini   = zeros([N,NA])+0.5;


for j = 1:NA
    for i = 1:N
        while delta_n(i,j) > toler
            c_n(i,j)     = (Pa/s_vo(i,j)).^n_ini(i,j);
            Qtn(i,j)     = ((qt(i)-svo(i,1))/Pa)*c_n(i,j);
            Fr(i,j)      = 100*1000*fs(i)/(qt(i)-svo(i,1));
            ic(i,j)      = ((3.47-log10(Qtn(i,j))).^2+(1.22+log10(Fr(i,j))).^2).^0.5;
            n_exp(i,j)   = min(1.0,0.381*(ic(i,j))+0.05*(s_vo(i,j)/Pa)-0.15);
            delta_n(i,j) = abs(n_exp(i,j) - n_ini(i,j));
            n_ini(i,j)   = n_exp(i,j);
        end
    end
end

%% Idriss and Boulanger Parameters
% Idriss and Boulanger (2014)
app_fc      = zeros([N,1]);
qcn         = (qt/Pa);
qc1n        = zeros([N,NA]);
d_qc1n      = zeros([N,NA]);
qc1ncs      = zeros([N,NA]);
cn          = zeros([N,NA]);
cn_g        = zeros([N,NA])+0.5;     % Initial guess
delta_cn    = zeros([N,NA])+0.5;     % Initial guess
toler_cn    = 0.01;

for j=1:NA
    for i = 1:N
        app_fc(i)   = 80*(ic(i)+cfc)-137;
        if app_fc(i) <0
            app_fc(i) = 0;
        elseif isnan(app_fc(i))
            app_fc(i) = NaN;
        elseif app_fc(i) >100
            app_fc(i)=100;
        end
        
        while delta_cn(i,j) > toler_cn
            qc1n(i,j)     = qcn(i)*cn_g(i,j);
            d_qc1n(i,j)   = (11.9+qc1n(i,j)/14.6)*exp(1.63-9.7/(app_fc(i)+2)-...
                (15.7/(2+app_fc(i))).^2);
            qc1ncs(i,j)   = qc1n(i,j) + d_qc1n(i,j);
            cn(i,j)       = min(1.7,(Pa/s_vo(i,j)).^(1.338-0.249*(qc1ncs(i,j)).^0.264));
            delta_cn(i,j) = abs(cn(i,j) - cn_g(i,j));
            cn_g(i,j)     = cn(i,j);
        end
    end
end

%% Evaluation of Liquefaction (B%I 14)
% Demand
CSR     = zeros([N,NA]);
rd      = zeros([N,NA]);
alpha_z = zeros([N,1]);
beta_z  = zeros([N,1]);

for i = 1:N
    alpha_z(i) = -1.012-1.126*sin(d(i)/11.73+5.133);
    beta_z(i)  = 0.106+0.118*sin(d(i)/11.28+5.142);
end

for j= 1:NA
    for i = 1:N
        rd(i,j)      = min(1.0,exp(alpha_z(i)+beta_z(i)*Mw(j,1)));
        CSR(i,j)     = 0.65*PGA(j,1)*rd(i,j)*svo(i,1)/s_vo(i,j);
    end
end

% Resistance
CRR_bl  = zeros([N,NA]);         %CRR for Mw = 7.5 and S'vo = 1atm
Csigma  = zeros([N,NA]);
Ksigma  = zeros([N,NA]);
MSF     = zeros([N,NA]);
CRR     = zeros([N,NA]);

if p_l==0.15
    C_o=2.8;
elseif p_l==0.5
    C_o=2.6;
end

for j=1:NA
    for i=1:N
        if ic(i,j)> 2.6 || d(i)<wt(j,1)
            CRR_bl(i,j) = NaN;
            MSF(i,j)    = NaN;
            Csigma(i,j) = NaN;
            Ksigma(i,j) = NaN;
            CRR(i,j)    = NaN;
        else
            CRR_bl(i,j) = exp(qc1ncs(i,j)/113+(qc1ncs(i,j)/1000)^2-(qc1ncs(i,j)/140)^3+((qc1ncs(i,j)/137)^4)-C_o);
            MSF(i,j)    = 1+((min(2.2,1.09+(qc1ncs(i,j)/180)^3))-1)*(8.64*exp(-Mw(j,1)/4)-1.325);
            Csigma(i,j) = min(0.3,max(0,1/(37.3-8.27*(qc1ncs(i,j)).^0.264)));
            Ksigma(i,j) = min(1.1,1-Csigma(i,j)*log(s_vo(i,j)/Pa));
            CRR(i,j)    = CRR_bl(i,j)*MSF(i,j)*Ksigma(i,j);
        end
    end
end

% Factor of Safety

FS = zeros([N,NA]);

for j=1:NA
    for i=1:N
        FS(i,j) = min(4,CRR(i,j)/CSR(i,j));
    end
end

% Probability of Liquefaction
CSR_bl      = zeros([N,NA]);         % CSR for Mw = 7.5 and S'vo = 1atm
Pl          = zeros([N,NA]);         % Probability of Liquefaction
Co          = 2.6;
Slnr        = 0.2;

for j=1:NA
    for i=1:N
        CSR_bl(i,j)  = CSR(i,j)./Ksigma(i,j);
        CSR_bl(i,j)  = CSR_bl(i,j)./MSF(i,j);
        if ic(i,j)>2.6 || d(i)<wt(j,1)
            Pl(i,j) = NaN;
        else
            Pl(i,j) =  100* (1-normcdf((qc1ncs(i,j)/113+(qc1ncs(i,j)/1000).^2-...
                (qc1ncs(i,j)/140).^3+(qc1ncs(i,j)/137).^4-Co-log(CSR_bl(i,j)))/Slnr,0,1));
        end
    end
end

%% Relative Denisty

% Correlation 1: Idriss and Boulanger (2008)
% Correlation 2: Jamiolkowsky (2001)
% Correlation 3: Kulhawy & Mayne

Dr_ib   = zeros([N,NA]);
Dr_j    = zeros([N,NA]);
%Qcn     = zeros([N,NA]);
Dr_km   = zeros([N,NA]);
bx      = 0.675;            %Compressibility factor (0.525,0.675,0.825
%for high, medium and low compressibility)
Qc      = 1.0   ;            % Compresibility factor (1.1 for high...
% 0.9 for low compressibility)
OCR    = 1.1 ;               % assumed OCR
Qocr   = OCR^0.18   ;              %OCR Factor (OCR^.18)
Qa     = 1.0;                % Age Factor 1.2+0.05log(t/100)

for j=1:NA
    for i = 1:N
        if ic(i,j)> 2.6 || d(i)<wt(j,1)
            Dr_ib(i,j) = NaN;
            Dr_j(i,j)  = NaN;
            Dr_km(i,j) = NaN;
        else
            Dr_ib(i,j)    = 100*(0.478*(qc1ncs(i,j)).^0.264 - 1.063);
            Dr_j(i,j)     = 100*(0.268*log((qt(i)/Pa)/(s_vo(i,j)/Pa).^0.5)-bx);
            Dr_km(i,j)    = 100*(((qt(i)/Pa)/(s_vo(i,j)/Pa).^0.5)/(305*Qc*Qocr*Qa))^0.5;
        end
    end
end

%% Evaluation of Liquefaction (Robertson 15)

% Demand
CSR_rob     = zeros([N,NA]);
rd_rob      = zeros([N,NA]);

for j= 1:NA
    for i = 1:N
        if d(i)<=9.15
            rd_rob(i,j)  = 1-0.00765*d(i);
        elseif d(i)<=23
            rd_rob(i,j)  = 1.174-0.0267*d(i);
        elseif d(i)<=30
            rd_rob(i,j)  = 0.744-0.008*d(i);
        else
            rd_rob(i,j)  = 0.5;
        end
        CSR_rob(i,j)     = 0.65*PGA(j,1)*rd_rob(i,j)*svo(i,1)/s_vo(i,j);
    end
end

% Resistance
% CRR_bl for Robertosn is given in Ic section
k_c         = zeros([N,NA]);
Qtncs       = zeros([N,NA]);
CRR_rob_bl  = zeros([N,NA]);
Ksigma_rob  = zeros([N,NA]);
CRR_rob     = zeros([N,NA]);
MSF_rob     = zeros([N,NA]);
factor_ks   = zeros([N,NA]);

for j= 1:NA
    for i=1:N
        if ic(i,j)<=2.50
            if ic(i,j) <= 1.64 || ic(i,j)>1.64 && ic(i,j)<2.36 && Fr(i,j)<0.5
                k_c(i,j) = 1.0;
            elseif ic(i,j)<=2.5
                k_c(i,j) = 5.58*ic(i,j).^3-0.403*ic(i,j).^4-21.63*ic(i,j).^2+33.75*ic(i,j)-17.88;
            end
        elseif ic(i,j)<=2.6
            k_c(i,j) = 6e-7*ic(i,j).^16.76;
        else
            k_c(i,j) =NaN;
        end
        Qtncs(i,j) = Qtn(i,j)*k_c(i,j);
        if Qtncs(i,j)<=50
            CRR_rob_bl(i,j) = 0.833*(Qtncs(i,j)/1000)+0.05;
        elseif Qtncs(i,j)<=160
            CRR_rob_bl(i,j) = 93*(Qtncs(i,j)/1000).^3+0.08;
        else
            CRR_rob_bl(i,j) = 0.5;
        end
    end
end


for j=1:NA
    for i=1:N
        if ic(i,j)> 2.6 || d(i)<wt(j,1)
            MSF_rob(i,j)    = NaN;
            Ksigma_rob(i,j) = NaN;
            factor_ks(i,j)  = NaN;
            CRR_rob(i,j)    = NaN;
            CRR_rob_bl(i,j) = NaN;
        else
            MSF_rob(i,j)    = 174/Mw(j,1).^2.56;
            if Dr_km(i,j)<=40
                factor_ks(i,j)  = 0.8;
            elseif Dr_km(i,j)<=80
                factor_ks(i,j)  =1-Dr_km(i,j)/200 ;
            else
                factor_ks(i,j)  = 0.6;
            end
            Ksigma_rob(i,j) = (s_vo(i,j)/Pa).^(factor_ks(i,j)-1);
            CRR_rob(i,j)    = CRR_rob_bl(i,j)*MSF_rob(i,j)*Ksigma_rob(i,j);
        end
    end
end

% Factor of Safety & Probability of Liquefaction

FS_rob = zeros([N,NA]);
Pl_rob          = zeros([N,NA]);         % Probability of Liquefaction

for j=1:NA
    for i=1:N
        FS_rob(i,j) = min(4,CRR_rob(i,j)/CSR_rob(i,j));
        Pl_rob(i,j) = 100/(1+(FS_rob(i,j)/0.9).^6.3);
    end
end

%% Shear Wave Velocity (based on McGann et al. (2014) for Christchurch)
Vs     = zeros([N,1]);
Vs_rob = zeros([N,NA]);

for j=1:NA
    for i = 1:N
        Vs(i,1) = 18.4* ((1000*qc(i,1))^0.144)* ((1000*fs(i,1))^0.0832) * (d(i,1)^0.278);
        Vs_rob(i,j) = (((10.^(0.55*ic(i,j)+1.68))*(qt(i)-svo(i)))/Pa)^0.5;
    end
end

%%%%%Added by JM, calculating Vs30
Vs30mc=zeros(NA,1);
Vs30rob=zeros(NA,1);
for ii=1:NA
    summc=0;
    sumrob=0;
    for jj= 1:N
        if isnan(Vs(jj,1))|| isnan(Vs_rob(jj,ii))
            % summc=summc;
            % sumrob=sumrob;
        else
            summc=summc+d_d/Vs(jj,1);
            sumrob=sumrob+d_d/Vs_rob(jj,ii);
        end
    end
    
    if d(N,1)<22
        summc=summc+(22-d(N,1))/Vs(N,1)+8/400;
        sumrob=sumrob+(22-d(N,1))/Vs_rob(N,ii)+8/400;
    else
        summc=summc+(30-d(N,1))/400;
        sumrob=sumrob+(30-d(N,1))/400;
    end
    
    Vs30mc(ii,1)=30/summc;
    Vs30rob(ii,1)=30/sumrob;
    
end

%%%%%%%%%%%%%%%%end of added by JM

%% Post liquefaction volumetric strain (IB08)
% BI14-IB08
gamma_lim_ib_bi = zeros([N,NA]);
F_alpha_ib_bi = zeros([N,NA]);
gamma_max_ib_bi = zeros([N,NA]);
vol_st_ib_bi = zeros([N,NA]);
post_liq_d_ib_bi = zeros([N,NA]);

for i=1:N
    for j=1:NA
        if isnan(FS(i,j)) || isnan(qc1ncs(i))
            gamma_lim_ib_bi(i,j) = NaN;
            F_alpha_ib_bi(i,j)   = NaN;
            gamma_max_ib_bi(i,j) = NaN;
            vol_st_ib_bi(i,j) = NaN;
        else
            gamma_lim_ib_bi(i,j) = max(0,1.859*(2.163-0.478*(qc1ncs(i,j)).^0.264).^3);
            F_alpha_ib_bi(i,j)   = -11.74+8.34*(min(168,max(69,qc1ncs(i,j)))).^0.264-1.371*(min(168,max(69,qc1ncs(i,j)))).^0.528;
            if FS(i,j)>=2
                gamma_max_ib_bi(i,j) = 0;
            elseif FS(i,j)>F_alpha_ib_bi(i,1)
                gamma_max_ib_bi(i,j) = min(gamma_lim_ib_bi(i,j),0.035*(2-FS(i,j))*(1-F_alpha_ib_bi(i,1))/(FS(i,j)-F_alpha_ib_bi(i,1)));
            else
                gamma_max_ib_bi(i,j) = gamma_lim_ib_bi(i,1);
            end
            vol_st_ib_bi(i,j) = 1.8*exp(2.551-1.147*(max(21,qc1ncs(i,j))).^0.264)*min(0.08,gamma_max_ib_bi(i,j));
        end
        
    end
end
for i=1:N
    for j=1:NA
        if isnan(FS(i,j)) || isnan(qc1ncs(i,j))
            post_liq_d_ib_bi(i,j) =  NaN;
        else
            if i==1
                post_liq_d_ib_bi(N-(i-1),j) = 100*vol_st_ib_bi(N-(i-1),j)*d_d;
            else
                post_liq_d_ib_bi(N-(i-1),j) = post_liq_d_ib_bi(N-i+2,j)+100*(vol_st_ib_bi(N-(i-1),j)*d_d);
            end
        end
    end
end
% R15-IB08
gamma_lim_rw_bi = zeros([N,NA]);
F_alpha_rw_bi = zeros([N,NA]);
gamma_max_rw_bi = zeros([N,NA]);
vol_st_rw_bi = zeros([N,NA]);
post_liq_d_rw_bi = zeros([N,NA]);

for i=1:N
    for j=1:NA
        if isnan(FS_rob(i,j)) || isnan(Qtncs(i))
            gamma_lim_rw_bi(i,j) = NaN;
            F_alpha_rw_bi(i,j)   = NaN;
            gamma_max_rw_bi(i,j) = NaN;
            vol_st_rw_bi(i,j) = NaN;
        else
            gamma_lim_rw_bi(i,j) = max(0,1.859*(2.163-0.478*(Qtncs(i,j)).^0.264).^3);
            F_alpha_rw_bi(i,j)   = -11.74+8.34*(min(168,max(69,Qtncs(i,j)))).^0.264-1.371*(min(168,max(69,Qtncs(i,j)))).^0.528;
            if FS_rob(i,j)>=2
                gamma_max_rw_bi(i,j) = 0;
            elseif FS_rob(i,j)>F_alpha_rw_bi(i,1)
                gamma_max_rw_bi(i,j) = min(gamma_lim_rw_bi(i,j),0.035*(2-FS_rob(i,j))*(1-F_alpha_rw_bi(i,1))/(FS_rob(i,j)-F_alpha_rw_bi(i,1)));
            else
                gamma_max_rw_bi(i,j) = gamma_lim_rw_bi(i,1);
            end
            vol_st_rw_bi(i,j) = 1.8*exp(2.551-1.147*(max(21,Qtncs(i,j))).^0.264)*min(0.08,gamma_max_rw_bi(i,j));
        end
        
    end
end
for i=1:N
    for j=1:NA
        if isnan(FS_rob(i,j)) || isnan(qc1ncs(i,j))
            post_liq_d_rw_bi(i,j) =  NaN;
        else
            if i==1
                post_liq_d_rw_bi(N-(i-1),j) = 100*vol_st_rw_bi(N-(i-1),j)*d_d;
            else
                post_liq_d_rw_bi(N-(i-1),j) = post_liq_d_rw_bi(N-i+2,j)+100*(vol_st_rw_bi(N-(i-1),j)*d_d);
            end
        end
    end
end

%% Post liquefaction volumetric strain (Z02)
vol_st_bi_z = zeros(N,NA);
%BI(14) - Z(02)
for i=1:N
    for j=1:NA
        if isnan(FS(i,j)) || isnan(qc1ncs(i,j))
            vol_st_rw_bi(i,j) = NaN;
        elseif FS(i,j)<=0.5
            vol_st_bi_z(i,j) = 102*qc1ncs(i,j).^-0.82;
        elseif FS(i,j)<=0.55
            vol_st_bi_z(i,j) = min(102*qc1ncs(i,j).^-0.82,( 102*qc1ncs(i,j).^-0.820)-(( 102*qc1ncs(i,j).^-0.820)-(2500*qc1ncs(i,j).^-1.440))*(FS(i,j)-0.50)/(0.05));
        elseif FS(i,j)<=0.60
            vol_st_bi_z(i,j) = min(102*qc1ncs(i,j).^-0.82,(2500*qc1ncs(i,j).^-1.440)-((2500*qc1ncs(i,j).^-1.440)-(2411*qc1ncs(i,j).^-1.450))*(FS(i,j)-0.55)/(0.05));
        elseif FS(i,j)<=0.65
            vol_st_bi_z(i,j) = min(102*qc1ncs(i,j).^-0.82,(2411*qc1ncs(i,j).^-1.450)-((2411*qc1ncs(i,j).^-1.450)-(2056*qc1ncs(i,j).^-1.435))*(FS(i,j)-0.60)/(0.05));
        elseif FS(i,j)<=0.70
            vol_st_bi_z(i,j) = min(102*qc1ncs(i,j).^-0.82,(2056*qc1ncs(i,j).^-1.435)-((2056*qc1ncs(i,j).^-1.435)-(1701*qc1ncs(i,j).^-1.420))*(FS(i,j)-0.65)/(0.05));
        elseif FS(i,j)<=0.75
            vol_st_bi_z(i,j) = min(102*qc1ncs(i,j).^-0.82,(1701*qc1ncs(i,j).^-1.420)-((1701*qc1ncs(i,j).^-1.420)-(1696*qc1ncs(i,j).^-1.440))*(FS(i,j)-0.70)/(0.05));
        elseif FS(i,j)<=0.80
            vol_st_bi_z(i,j) = min(102*qc1ncs(i,j).^-0.82,(1696*qc1ncs(i,j).^-1.440)-((1696*qc1ncs(i,j).^-1.440)-(1690*qc1ncs(i,j).^-1.460))*(FS(i,j)-0.75)/(0.05));
        elseif FS(i,j)<=0.85
            vol_st_bi_z(i,j) = min(102*qc1ncs(i,j).^-0.82,(1690*qc1ncs(i,j).^-1.460)-((1690*qc1ncs(i,j).^-1.460)-(1560*qc1ncs(i,j).^-1.470))*(FS(i,j)-0.80)/(0.05));
        elseif FS(i,j)<=0.90
            vol_st_bi_z(i,j) = min(102*qc1ncs(i,j).^-0.82,(1560*qc1ncs(i,j).^-1.470)-((1560*qc1ncs(i,j).^-1.470)-(1430*qc1ncs(i,j).^-1.480))*(FS(i,j)-0.85)/(0.05));
        elseif FS(i,j)<=0.95
            vol_st_bi_z(i,j) = min(102*qc1ncs(i,j).^-0.82,(1430*qc1ncs(i,j).^-1.480)-((1430*qc1ncs(i,j).^-1.480)-( 150*qc1ncs(i,j).^-1.070))*(FS(i,j)-0.90)/(0.05));
        elseif FS(i,j)<=1.00
            vol_st_bi_z(i,j) = min(102*qc1ncs(i,j).^-0.82,( 150*qc1ncs(i,j).^-1.070)-(( 150*qc1ncs(i,j).^-1.070)-(  64*qc1ncs(i,j).^-0.930))*(FS(i,j)-0.95)/(0.05));
        elseif FS(i,j)<=1.05
            vol_st_bi_z(i,j) = min(102*qc1ncs(i,j).^-0.82,(  64*qc1ncs(i,j).^-0.930)-((  64*qc1ncs(i,j).^-0.930)-(  31*qc1ncs(i,j).^-0.820))*(FS(i,j)-1.00)/(0.05));
        elseif FS(i,j)<=1.10
            vol_st_bi_z(i,j) = min(102*qc1ncs(i,j).^-0.82,(  31*qc1ncs(i,j).^-0.820)-((  31*qc1ncs(i,j).^-0.820)-(  11*qc1ncs(i,j).^-0.650))*(FS(i,j)-1.05)/(0.05));
        elseif FS(i,j)<=1.15
            vol_st_bi_z(i,j) = min(102*qc1ncs(i,j).^-0.82,(  11*qc1ncs(i,j).^-0.650)-((  11*qc1ncs(i,j).^-0.650)-(10.4*qc1ncs(i,j).^-0.670))*(FS(i,j)-1.10)/(0.05));
        elseif FS(i,j)<=1.20
            vol_st_bi_z(i,j) = min(102*qc1ncs(i,j).^-0.82,(10.4*qc1ncs(i,j).^-0.670)-((10.4*qc1ncs(i,j).^-0.670)-( 9.7*qc1ncs(i,j).^-0.690))*(FS(i,j)-1.15)/(0.05));
        elseif FS(i,j)<=1.25
            vol_st_bi_z(i,j) = min(102*qc1ncs(i,j).^-0.82,( 9.7*qc1ncs(i,j).^-0.690)-(( 9.7*qc1ncs(i,j).^-0.690)-( 8.7*qc1ncs(i,j).^-0.700))*(FS(i,j)-1.20)/(0.05));
        elseif FS(i,j)<=1.30
            vol_st_bi_z(i,j) = min(102*qc1ncs(i,j).^-0.82,( 8.7*qc1ncs(i,j).^-0.700)-(( 8.7*qc1ncs(i,j).^-0.700)-( 7.6*qc1ncs(i,j).^-0.710))*(FS(i,j)-1.25)/(0.05));
        elseif FS(i,j)<=1.40
            vol_st_bi_z(i,j) = min(102*qc1ncs(i,j).^-0.82,( 7.6*qc1ncs(i,j).^-0.710)-(( 7.6*qc1ncs(i,j).^-0.710)-( 6.0*qc1ncs(i,j).^-0.683))*(FS(i,j)-1.30)/(0.10));
        elseif FS(i,j)<=1.50
            vol_st_bi_z(i,j) = min(102*qc1ncs(i,j).^-0.82,( 6.0*qc1ncs(i,j).^-0.683)-(( 6.0*qc1ncs(i,j).^-0.683)-( 5.0*qc1ncs(i,j).^-0.683))*(FS(i,j)-1.40)/(0.10));
        elseif FS(i,j)<=1.60
            vol_st_bi_z(i,j) = min(102*qc1ncs(i,j).^-0.82,( 5.0*qc1ncs(i,j).^-0.683)-(( 5.0*qc1ncs(i,j).^-0.683)-( 4.0*qc1ncs(i,j).^-0.683))*(FS(i,j)-1.50)/(0.10));
        elseif FS(i,j)<=1.70
            vol_st_bi_z(i,j) = min(102*qc1ncs(i,j).^-0.82,( 4.0*qc1ncs(i,j).^-0.683)-(( 4.0*qc1ncs(i,j).^-0.683)-( 3.0*qc1ncs(i,j).^-0.683))*(FS(i,j)-1.60)/(0.10));
        elseif FS(i,j)<=1.80
            vol_st_bi_z(i,j) = min(102*qc1ncs(i,j).^-0.82,( 3.0*qc1ncs(i,j).^-0.683)-(( 3.0*qc1ncs(i,j).^-0.683)-( 2.0*qc1ncs(i,j).^-0.683))*(FS(i,j)-1.70)/(0.10));
        elseif FS(i,j)<=1.90
            vol_st_bi_z(i,j) = min(102*qc1ncs(i,j).^-0.82,( 2.0*qc1ncs(i,j).^-0.683)-(( 2.0*qc1ncs(i,j).^-0.683)-( 1.0*qc1ncs(i,j).^-0.683))*(FS(i,j)-1.80)/(0.10));
        elseif FS(i,j)<=2.0
            vol_st_bi_z(i,j) = min(102*qc1ncs(i,j).^-0.82,( 1.0*qc1ncs(i,j).^-0.683)-(( 1.0*qc1ncs(i,j).^-0.683)*(FS(i,j)-1.80)/(0.10)));
        else
            vol_st_bi_z(i,j) = 0.0;
        end
    end
end

post_liq_d_bi_z=zeros(N,NA);
for i=1:N
    for j=1:NA
        if isnan(FS(i,j)) || isnan(qc1ncs(i,j))
            %post_liq_d_bi_z(i,j) =  NaN;
            
            if i==1
                post_liq_d_bi_z(N-(i-1),j) = vol_st_bi_z(N-(i-1),j)*d_d;
            else
                post_liq_d_bi_z(N-(i-1),j) = post_liq_d_bi_z(N-i+2,j);
            end
            
        else
            if i==1
                post_liq_d_bi_z(N-(i-1),j) = vol_st_bi_z(N-(i-1),j)*d_d;
            else
                post_liq_d_bi_z(N-(i-1),j) = post_liq_d_bi_z(N-i+2,j)+(vol_st_bi_z(N-(i-1),j)*d_d);
            end
        end
    end
end
%%%Added for LSN Versus depth
post_liq_d_bi_z_LSN=zeros(N,NA);
for i=1:N
    for j=1:NA
        if isnan(FS(i,j)) || isnan(qc1ncs(i,j))
            %post_liq_d_bi_z_LSN(i,j) =  NaN;
            if i==1
                post_liq_d_bi_z_LSN(N-(i-1),j) = vol_st_bi_z(N-(i-1),j)*d_d*1000/d(N-(i-1))/100;
                %post_liq_d_bi_z_LSN(N-(i-1),j) = post_liq_d_bi_z_LSN(N-i+2,j);
            else
                post_liq_d_bi_z_LSN(N-(i-1),j) = post_liq_d_bi_z_LSN(N-i+2,j);
            end
        else
            if i==1
                post_liq_d_bi_z_LSN(N-(i-1),j) = vol_st_bi_z(N-(i-1),j)*d_d*1000/d(N-(i-1))/100;
            else
                post_liq_d_bi_z_LSN(N-(i-1),j) = post_liq_d_bi_z_LSN(N-i+2,j)+(vol_st_bi_z(N-(i-1),j)*d_d*1000/d(N-(i-1))/100);
            end
        end
    end
end




%%%%%End of Added for LSN Versus depth
%%%%%%%Adding calculation of LPI
LPI=zeros(size(Mw,1),1);
for i=1:NA
    sum=0;
    for j=1:N
        
        if isnan(FS(j,i))|| isnan(qc1ncs(j,i))
            %sum=sum;%post_liq_d_bi_z(i,j) =  NaN;
        else
            if FS(j,i)<1
                F1=1-FS(j,i);
            else
                F1=0;
            end
            sum=sum+F1*(10-0.5*d(j))*d_d;
        end
    end
    LPI(i,1)=sum;
end

%%%%%%%%%%%%end LPI

%%%%%%%Adding calculation of LSN
LSN=zeros(size(Mw,1),1);
LBM1=zeros(size(Mw,1),1);
LSNO2=zeros(size(Mw,1),1);
for i=1:NA
    sum1=0;
    sum2=0;
    sum3=0;
    for j=1:N
        
        if isnan(FS(j,i))|| isnan(qc1ncs(j,i))
            %sum1=sum1;%post_liq_d_bi_z(i,j) =  NaN;
            %sum2=sum2;
        else
            wz1=1;
            if Df(i,1)>d(j)
                wz1=0;
            end
            ev=volset(FS(j,i),qc1ncs(j,i));
            sum1=sum1+1000*ev/100/d(j)*d_d;
            DR_prom=1/3*(Dr_ib(j,i)+Dr_j(j,i)+Dr_km(j,i));
            sh=shearset(FS(j,i),DR_prom);
            sum2=sum2+1000*sh/100/d(j)*wz1*d_d;
            sum3=sum3+1000*vol_st_bi_z(j,i)/100/d(j)*d_d;
        end
    end
    LSN(i,1)=sum1;
    LBM1(i,1)=sum2;
    LSNO2(i,1)=sum3;
end

%%%%%%%Adding calculation of LBM6
LBM1N=zeros(size(Mw,1),1);

for i=1:NA
    sum2=0;
    for j=1:N
        
        if isnan(FS(j,i))|| isnan(qc1ncs(j,i))
        else
            wz1=1;
            if Df(i,1)>d(j)
                wz1=0;
            end
            
            DR_prom=1/3*(Dr_ib(j,i)+Dr_j(j,i)+Dr_km(j,i));
            sh=shearset(FS(j,i),DR_prom);
            sum2=sum2+100*sh/100/d(j)*wz1*d_d;
            
        end
    end
    
    LBM1N(i,1)=sum2;
    
end
LBS = LBM1N;
%%%%%%%%%%%%end LBM1N

%%%%%%%%%%%%%% HL (B&I, robert)
idx = find(FS<1);
HL_BI = size(idx,1) * d_d;

idx = find(FS_rob<1);
HL_R15 = size(idx,1) * d_d;



vol_st_rw_z=zeros(N,NA);
%R(15) - Z(02)
for i=1:N
    for j=1:NA
        if isnan(FS_rob(i,j)) || isnan(Qtncs(i))
            vol_st_rw_z(i,j) = NaN;
        elseif FS_rob(i,j)<=0.5
            vol_st_rw_z(i,j) = 102*Qtncs(i,j).^-0.82;
        elseif FS_rob(i,j)<=0.55
            vol_st_rw_z(i,j) = min(102*Qtncs(i,j).^-0.82,( 102*Qtncs(i,j).^-0.820)-(( 102*Qtncs(i,j).^-0.820)-(2500*Qtncs(i,j).^-1.440))*(FS_rob(i,j)-0.50)/(0.05));
        elseif FS_rob(i,j)<=0.60
            vol_st_rw_z(i,j) = min(102*Qtncs(i,j).^-0.82,(2500*Qtncs(i,j).^-1.440)-((2500*Qtncs(i,j).^-1.440)-(2411*Qtncs(i,j).^-1.450))*(FS_rob(i,j)-0.55)/(0.05));
        elseif FS_rob(i,j)<=0.65
            vol_st_rw_z(i,j) = min(102*Qtncs(i,j).^-0.82,(2411*Qtncs(i,j).^-1.450)-((2411*Qtncs(i,j).^-1.450)-(2056*Qtncs(i,j).^-1.435))*(FS_rob(i,j)-0.60)/(0.05));
        elseif FS_rob(i,j)<=0.70
            vol_st_rw_z(i,j) = min(102*Qtncs(i,j).^-0.82,(2056*Qtncs(i,j).^-1.435)-((2056*Qtncs(i,j).^-1.435)-(1701*Qtncs(i,j).^-1.420))*(FS_rob(i,j)-0.65)/(0.05));
        elseif FS_rob(i,j)<=0.75
            vol_st_rw_z(i,j) = min(102*Qtncs(i,j).^-0.82,(1701*Qtncs(i,j).^-1.420)-((1701*Qtncs(i,j).^-1.420)-(1696*Qtncs(i,j).^-1.440))*(FS_rob(i,j)-0.70)/(0.05));
        elseif FS_rob(i,j)<=0.80
            vol_st_rw_z(i,j) = min(102*Qtncs(i,j).^-0.82,(1696*Qtncs(i,j).^-1.440)-((1696*Qtncs(i,j).^-1.440)-(1690*Qtncs(i,j).^-1.460))*(FS_rob(i,j)-0.75)/(0.05));
        elseif FS_rob(i,j)<=0.85
            vol_st_rw_z(i,j) = min(102*Qtncs(i,j).^-0.82,(1690*Qtncs(i,j).^-1.460)-((1690*Qtncs(i,j).^-1.460)-(1560*Qtncs(i,j).^-1.470))*(FS_rob(i,j)-0.80)/(0.05));
        elseif FS_rob(i,j)<=0.90
            vol_st_rw_z(i,j) = min(102*Qtncs(i,j).^-0.82,(1560*Qtncs(i,j).^-1.470)-((1560*Qtncs(i,j).^-1.470)-(1430*Qtncs(i,j).^-1.480))*(FS_rob(i,j)-0.85)/(0.05));
        elseif FS_rob(i,j)<=0.95
            vol_st_rw_z(i,j) = min(102*Qtncs(i,j).^-0.82,(1430*Qtncs(i,j).^-1.480)-((1430*Qtncs(i,j).^-1.480)-( 150*Qtncs(i,j).^-1.070))*(FS_rob(i,j)-0.90)/(0.05));
        elseif FS_rob(i,j)<=1.00
            vol_st_rw_z(i,j) = min(102*Qtncs(i,j).^-0.82,( 150*Qtncs(i,j).^-1.070)-(( 150*Qtncs(i,j).^-1.070)-(  64*Qtncs(i,j).^-0.930))*(FS_rob(i,j)-0.95)/(0.05));
        elseif FS_rob(i,j)<=1.05
            vol_st_rw_z(i,j) = min(102*Qtncs(i,j).^-0.82,(  64*Qtncs(i,j).^-0.930)-((  64*Qtncs(i,j).^-0.930)-(  31*Qtncs(i,j).^-0.820))*(FS_rob(i,j)-1.00)/(0.05));
        elseif FS_rob(i,j)<=1.10
            vol_st_rw_z(i,j) = min(102*Qtncs(i,j).^-0.82,(  31*Qtncs(i,j).^-0.820)-((  31*Qtncs(i,j).^-0.820)-(  11*Qtncs(i,j).^-0.650))*(FS_rob(i,j)-1.05)/(0.05));
        elseif FS_rob(i,j)<=1.15
            vol_st_rw_z(i,j) = min(102*Qtncs(i,j).^-0.82,(  11*Qtncs(i,j).^-0.650)-((  11*Qtncs(i,j).^-0.650)-(10.4*Qtncs(i,j).^-0.670))*(FS_rob(i,j)-1.10)/(0.05));
        elseif FS_rob(i,j)<=1.20
            vol_st_rw_z(i,j) = min(102*Qtncs(i,j).^-0.82,(10.4*Qtncs(i,j).^-0.670)-((10.4*Qtncs(i,j).^-0.670)-( 9.7*Qtncs(i,j).^-0.690))*(FS_rob(i,j)-1.15)/(0.05));
        elseif FS_rob(i,j)<=1.25
            vol_st_rw_z(i,j) = min(102*Qtncs(i,j).^-0.82,( 9.7*Qtncs(i,j).^-0.690)-(( 9.7*Qtncs(i,j).^-0.690)-( 8.7*Qtncs(i,j).^-0.700))*(FS_rob(i,j)-1.20)/(0.05));
        elseif FS_rob(i,j)<=1.30
            vol_st_rw_z(i,j) = min(102*Qtncs(i,j).^-0.82,( 8.7*Qtncs(i,j).^-0.700)-(( 8.7*Qtncs(i,j).^-0.700)-( 7.6*Qtncs(i,j).^-0.710))*(FS_rob(i,j)-1.25)/(0.05));
        elseif FS_rob(i,j)<=1.40
            vol_st_rw_z(i,j) = min(102*Qtncs(i,j).^-0.82,( 7.6*Qtncs(i,j).^-0.710)-(( 7.6*Qtncs(i,j).^-0.710)-( 6.0*Qtncs(i,j).^-0.683))*(FS_rob(i,j)-1.30)/(0.10));
        elseif FS_rob(i,j)<=1.50
            vol_st_rw_z(i,j) = min(102*Qtncs(i,j).^-0.82,( 6.0*Qtncs(i,j).^-0.683)-(( 6.0*Qtncs(i,j).^-0.683)-( 5.0*Qtncs(i,j).^-0.683))*(FS_rob(i,j)-1.40)/(0.10));
        elseif FS_rob(i,j)<=1.60
            vol_st_rw_z(i,j) = min(102*Qtncs(i,j).^-0.82,( 5.0*Qtncs(i,j).^-0.683)-(( 5.0*Qtncs(i,j).^-0.683)-( 4.0*Qtncs(i,j).^-0.683))*(FS_rob(i,j)-1.50)/(0.10));
        elseif FS_rob(i,j)<=1.70
            vol_st_rw_z(i,j) = min(102*Qtncs(i,j).^-0.82,( 4.0*Qtncs(i,j).^-0.683)-(( 4.0*Qtncs(i,j).^-0.683)-( 3.0*Qtncs(i,j).^-0.683))*(FS_rob(i,j)-1.60)/(0.10));
        elseif FS_rob(i,j)<=1.80
            vol_st_rw_z(i,j) = min(102*Qtncs(i,j).^-0.82,( 3.0*Qtncs(i,j).^-0.683)-(( 3.0*Qtncs(i,j).^-0.683)-( 2.0*Qtncs(i,j).^-0.683))*(FS_rob(i,j)-1.70)/(0.10));
        elseif FS_rob(i,j)<=1.90
            vol_st_rw_z(i,j) = min(102*Qtncs(i,j).^-0.82,( 2.0*Qtncs(i,j).^-0.683)-(( 2.0*Qtncs(i,j).^-0.683)-( 1.0*Qtncs(i,j).^-0.683))*(FS_rob(i,j)-1.80)/(0.10));
        elseif FS_rob(i,j)<=2.0
            vol_st_rw_z(i,j) = min(102*Qtncs(i,j).^-0.82,( 1.0*Qtncs(i,j).^-0.683)-(( 1.0*Qtncs(i,j).^-0.683)*(FS_rob(i,j)-1.80)/(0.10)));
        else
            vol_st_rw_z(i,j) = 0.0;
        end
    end
end


post_liq_d_rw_z = zeros(N,NA);
for i=1:N
    for j=1:NA
        if isnan(FS_rob(i,j)) || isnan(Qtncs(i,j))
            post_liq_d_rw_z(i,j) =  NaN;
        else
            if i==1
                post_liq_d_rw_z(N-(i-1),j) = vol_st_rw_z(N-(i-1),j)*d_d;
            else
                post_liq_d_rw_z(N-(i-1),j) = post_liq_d_rw_z(N-i+2,j)+(vol_st_rw_z(N-(i-1),j)*d_d);
            end
        end
    end
end


function [ev] = volset(FS,qc1ncs)

ev=0;
if FS<=0.5 && qc1ncs>=33 && qc1ncs<=200
    ev=102*(qc1ncs)^-0.82;
end

if FS>0.5 && FS<0.65 && qc1ncs>=33 && qc1ncs<=147
    ev=102*(qc1ncs)^-0.82;
end

if FS>0.5 && FS<0.65 && qc1ncs>147 && qc1ncs<=200
    ev=2411*(qc1ncs)^-1.45;
end

if FS>=0.65 && FS<0.75 && qc1ncs>=33 && qc1ncs<=110
    ev=102*(qc1ncs)^-0.82;
end

if FS>=0.65 && FS<0.75 && qc1ncs>110 && qc1ncs<=200
    ev=1701*(qc1ncs)^-1.42;
end

if FS>=0.75 && FS<0.85 && qc1ncs>=33 && qc1ncs<=80
    ev=102*(qc1ncs)^-0.82;
end

if FS>=0.75 && FS<0.85 && qc1ncs>80 && qc1ncs<=200
    ev=1690*(qc1ncs)^-1.46;
end

if FS>=0.85 && FS<0.95 && qc1ncs>=33 && qc1ncs<=60
    ev=102*(qc1ncs)^-0.82;
end
if FS>=0.85 && FS<0.95 && qc1ncs>60 && qc1ncs<=200
    ev=1430*(qc1ncs)^-1.48;
end

if FS>=0.95 && FS<1.05 && qc1ncs>=33 && qc1ncs<=200
    ev=64*(qc1ncs)^-0.93;
end

if FS>=1.05 && FS<1.15 && qc1ncs>=33 && qc1ncs<=200
    ev=11*(qc1ncs)^-0.65;
end

if FS>=1.15 && FS<1.25 && qc1ncs>=33 && qc1ncs<=200
    ev=9.7*(qc1ncs)^-0.69;
end

if FS>=1.25 && FS<1.35 && qc1ncs>=33 && qc1ncs<=200
    ev=7.6*(qc1ncs)^-0.71;
end

if FS>=1.35 && FS<2.0 && qc1ncs>=33 && qc1ncs<=200
    ev=7.6*(qc1ncs)^-0.71;
end

if FS>=2
    ev=0;
    
end


function [sh] = shearset(FS,DR)

sh=0;
if DR<=45 && FS<0.81
    sh=51.2;
end

if DR<=45 && FS>=0.81 && FS<=1.0
    sh=250*(1-FS)+3.5;
end

if DR<=45 && FS>=1.0 && FS<=2.0
    sh=3.31*FS^-7.97;
end


if DR>45&& DR<=55 && FS<=0.72
    sh=34.1;
end

if DR>45&& DR<=55 && FS>=0.72 && FS<=2.0
    sh=4.22*FS^-6.39;
end

if DR>55&& DR<=65 && FS<=0.66
    sh=22.7;
end

if DR>55&& DR<=65 && FS>=0.66 && FS<=2.0
    sh=3.58*FS^-4.42;
end

if DR>65&& DR<=75 && FS<=0.59
    sh=14.5;
end

if DR>65&& DR<=75 && FS>=0.59 && FS<=2.0
    sh=3.20*FS^-2.89;
end

if DR>75&& DR<=85 && FS<=0.56
    sh=10;
end

if DR>75&& DR<=85 && FS>=0.56 && FS<=2.0
    sh=3.22*FS^-2.08;
end

if DR>85&& DR<=95 && FS<=0.70
    sh=6.2;
end

if DR>85&& DR<=95 && FS>=0.70 && FS<=2.0
    sh=3.26*FS^-1.80;
end







