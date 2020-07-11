function [lny,sigma,tau,phi]=Kuehn2020(To,M,Rrup,Ztor,Vs30,mechanism,region,alfaBackArc,alfaNankai,Z10,Z25,Nepist)

st          = dbstack;
isadmisible = isIMadmisible(To,st(1).name,[0 10],[nan nan],[nan nan],[nan nan]);
if isadmisible==0
    lny   = nan(size(M));
    sigma = nan(size(M));
    tau   = nan(size(M));
    phi   = nan(size(M));
    return
end


Z10 = Z10*1000; %Z1.0 must be input in km
Z25 = Z25*1000; %Z2.5 must be input in km 

% mechanism flag and Magnitude breaking point
switch mechanism
    case 'interface'
        Fs=0;
        switch region
            case 'global',               MB=7.8;   regID='All';
            case 'aleutian',             MB=8.0;   regID='AK';
            case 'alaska',               MB=8.0;   regID='AK';
            case 'cascadia',             MB=7.56;  regID='CASC';
            case 'central_america_s',    MB=7.5;   regID='CAM';
            case 'central_america_n',    MB=7.45;  regID='CAM';
            case 'japan_pac',            MB=8.31;  regID='JP';
            case 'japan_phi',            MB=7.28;  regID='JP';
            case 'new_zealand_n',        MB=7.8;   regID='NZ';
            case 'new_zealand_s',        MB=7.8;   regID='NZ';
            case 'south_america_n',      MB=8.45;  regID='SA';
            case 'south_america_s',      MB=8.45;  regID='SA';
            case 'taiwan_w',             MB=8.0;   regID='TW';
            case 'taiwan_e',             MB=8.0;   regID='TW';
        end
    case 'intraslab'
        Fs=1;
        switch region
            case 'global',               MB=7.6;   regID='All';
            case 'aleutian',             MB=7.98;  regID='AK';
            case 'alaska',               MB=7.2;   regID='AK';
            case 'cascadia',             MB=7.2;   regID='CASC';
            case 'central_america_s',    MB=7.6;   regID='CAM';
            case 'central_america_n',    MB=7.4;   regID='CAM';
            case 'japan_pac',            MB=7.65;  regID='JP';
            case 'japan_phi',            MB=7.55;  regID='JP';
            case 'new_zealand_n',        MB=7.6;   regID='NZ';
            case 'new_zealand_s',        MB=7.4;   regID='NZ';
            case 'south_america_n',      MB=7.3;   regID='SA';
            case 'south_america_s',      MB=7.25;  regID='SA';
            case 'taiwan_w',             MB=7.7;   regID='TW';
            case 'taiwan_e',             MB=7.7;   regID='TW';
        end
end

Nscen = length(Rrup);
zer = zeros(Nscen,1);
switch regID
    case 'All',  R1 = alfaBackArc * Rrup;  R3 = zer;
    case 'AK',   R1 = zer;                 R3 = zer;
    case 'CASC', R1 = zer;                 R3 = zer;
    case 'CAM',  R1 = alfaBackArc * Rrup;  R3 = zer;
    case 'JP',   R1 = alfaBackArc * Rrup;  R3 = alfaNankai* Rrup;
    case 'NZ',   R1 = zer;                 R3 = zer;
    case 'SA',   R1 = alfaBackArc * Rrup;  R3 = zer;
    case 'TW',   R1 = zer;                 R3 = zer;
end
R2 = Rrup - R1 - R3;
Fx = ones(size(Rrup));
Fx(and(R1==0,R2==0))=0;
Fx(and(R2==0,R3==0))=0;
Fx(and(R1==0,R3==0))=0;

if To>=0
    To      = max(To,0.001);  % PGA is associated to To=0.01;
end
period  = [-1 0.001 0.01,0.02,0.03,0.05,0.075,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.75,1,1.5,2,3,4,5,7.5,10];
T_lo    = max(period(period<=To));
T_hi    = min(period(period>=To));
index   = find(abs((period - T_lo)) < 1e-6);    % Identify the period

Cpga    = Kuehn2020_coeff(2,regID);
PGA1100 = exp(Kuehn2020_PGA(Cpga,M,MB,Rrup,1100,Ztor,Fs,regID,Fx,R1,R2,R3,Z10,Z25,-999));

if T_lo==T_hi
    C             = Kuehn2020_coeff(index,regID);
    [lny,tau,phi] = gmpe(index,C,Cpga,M,MB,Rrup,Vs30,Ztor,Fs,regID,Fx,R1,R2,R3,Z10,Z25,PGA1100);
else
    C             = Kuehn2020_coeff(index  ,regID);
    Cp1           = Kuehn2020_coeff(index+1,regID);
    [lny_lo,tau_lo,phi_lo] = gmpe(index,    C,Cpga,M,MB,Rrup,Vs30,Ztor,Fs,regID,Fx,R1,R2,R3,Z10,Z25,PGA1100);
    [lny_hi,tau_hi,phi_hi] = gmpe(index+1,Cp1,Cpga,M,MB,Rrup,Vs30,Ztor,Fs,regID,Fx,R1,R2,R3,Z10,Z25,PGA1100);
    x          = log([T_lo;T_hi]);
    Y_sa       = [lny_lo,lny_hi]';
    Y_sigma2   = [tau_lo,tau_hi]';
    Y_sigma1   = [phi_lo,phi_hi]';
    lny        = interp1(x,Y_sa,log(To))';
    phi        = interp1(x,Y_sigma2,log(To))';
    tau        = interp1(x,Y_sigma1,log(To))';
end

sigma_aleatory = sqrt(tau.^2+phi.^2); 

if Nepist==0
    sig_epistemic=zer;
else
    Tstring={'T-1.00','T00.00','T00.01','T00.02','T00.03','T00.05','T00.07','T00.10','T00.15','T00.20','T00.25','T00.30','T00.40','T00.50','T00.75','T01.00','T01.50','T02.00','T03.00','T04.00','T05.00','T07.50','T10.00'};
    fname = sprintf('posterior_coefficients_KBCG20_%s.csv',Tstring{index});
    data  = csvread(fname,1,1,[1 1 Nepist 179]);
    med   = zeros(Nscen,Nepist);
    for i=1:Nepist
        C        = Kuehn2020_coeff(index,regID,data(i,:));
        Cpga     = Kuehn2020_coeff(2    ,regID,data(i,:));
        med(:,i) = gmpe(index,C,Cpga,M,MB,Rrup,Vs30,Ztor,Fs,regID,Fx,R1,R2,R3,Z10,Z25,PGA1100);
    end
    sig_epistemic = std(med,0,2);
end

sigma  = sqrt(sigma_aleatory.^2+sig_epistemic.^2);

function[lny,tau,phi]=gmpe(index,C,Cpga,M,MB,Rrup,Vs30,Ztor,Fs,regID,Fx,R1,R2,R3,Z10,Z25,PGA1100)
period = [-1 0.001 0.01,0.02,0.03,0.05,0.075,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.75,1,1.5,2,3,4,5,7.5,10];
T      = period(index);
dzb_if =C(1 );
dzb_sl =C(2 );
nft1   =C(3 );
nft2   =C(4 );
phi    =C(5 );
tau    =C(6 );
th1_if =C(7 );
th1_sl =C(8 );
th2_if =C(9 );
th2_sl =C(10);
th3    =C(11);
th4_if =C(12);
th4_sl =C(13);
th5    =C(14);
th6_1  =C(15);
th6_2  =C(16);
th6_3  =C(17);
th6_x1 =C(18);
th6_x2 =C(19);
th6_x3 =C(20);
th6_xc =C(21);
th7    =C(22);
th9_if =C(23);
th9_sl =C(24);
k1     =C(25);
k2     =C(26);
thz1   =C(27);
thz2   =C(28);
thz3   =C(29);
thz4   =C(30);
th11   =C(31);
th12   =C(32);

% Magnitude scaling term
deltaM   = 0.1;
deltaM_B = 0;
Mref     = 6;
deltaZ   = 1;
th10     = 0;
ZB_if    = 30;
ZB_sl    = 80;

if Fs==0
    if and(T>1,T<3)
        deltaM_B=interp1(log([1 3]),[0 -0.4],log(T),'linear');
    elseif T>=3
        deltaM_B=-0.4;
    end
    MB = MB + deltaM_B;
end

% MAGNITUDE TERM
switch Fs
    case 0, fmag = lh(M,MB,th4_if*(MB - Mref),th4_if,th5,deltaM);
    case 1, fmag = lh(M,MB,th4_sl*(MB - Mref),th4_sl,th5,deltaM);
end

% GEOMETRICAL ATTENUATION TERM
hM = 10.^(nft1+nft2*(M-6));
switch Fs
    case 0, fgeom = (th2_if + th3*M).*log(Rrup+hM);
    case 1, fgeom = (th2_sl + th3*M).*log(Rrup+hM);
end

% DEPTH TERM
Zif_ref = 15;
Zsl_ref = 50;
switch Fs
    case 0, fdepth = lh(Ztor,ZB_if + dzb_if,th9_if*(ZB_if + dzb_if-Zif_ref),th9_if,th10,deltaZ);
    case 1, fdepth = lh(Ztor,ZB_sl + dzb_sl,th9_sl*(ZB_sl + dzb_sl-Zsl_ref),th9_sl,th10,deltaZ);
end

%% SITE TERM
c	= 1.88;
n	= 1.18;
if Vs30<=k1
    fsite = th7*log(Vs30/k1) + k2*(log(PGA1100+c*(Vs30/k1)^n)-log(PGA1100+c));
else
    fsite = (th7+k2*n)*log(Vs30/k1);
end

%% ANAELASTIC ATTENUATION TERM
switch regID
    case 'All'
        fattn = Fx.*(th6_xc+th6_x1*R1+th6_x2*R2+th6_x3*R3)+(1-Fx).*(th6_1*R1+th6_2*R2+th6_3*R3);    
    case 'JP'
        fattn = Fx.*(th6_xc+th6_x1*R1+th6_x2*R2+th6_x3*R3)+(1-Fx).*(th6_1*R1+th6_2*R2+th6_3*R3);
    case {'CAM','SA'}
        fattn = Fx.*(th6_xc+th6_x1*R1+th6_x2*R2)+(1-Fx).*(th6_1*R1+th6_2*R2);
    case {'AK','CASC','NZ','TW'}
        fattn = th6_2*Rrup;
    otherwise
        fattn=0;
end

%% BASIN DEPTH TERM
aux    = (log(Vs30)-thz3)/thz4;
lnZref = thz1+(thz2-thz1)*exp(aux)/(1+exp(aux));

switch regID
    case 'CASC'
        delta_lnZ = log(Z25)-lnZref;
    case {'NZ','TW'}
        delta_lnZ = log(Z10)-lnZref;
    case 'JP'
        delta_lnZ = log(Z25)-lnZref;
    otherwise
        th11      = 0;
        th12      = 0;
        delta_lnZ = 0;
end

fbasin = th11 + th12 *delta_lnZ;

%% BASE COEFFICIENT
lny = (1-Fs)*th1_if + Fs*th1_sl + fmag +fgeom + fdepth + fattn + fsite + fbasin;

%% SPECTRAL SHAPE CORRECTION AT LOW PERIODS
if T<=0.1
    lnPGA   = Kuehn2020_PGA(Cpga,M,MB,Rrup,Vs30,Ztor,Fs,regID,Fx,R1,R2,R3,Z10,Z25,PGA1100);
    ind     = lny<lnPGA;
    lny(ind)= lnPGA(ind);
end

function[lnPGA]=Kuehn2020_PGA(Cpga,M,MB,Rrup,Vs30,Ztor,Fs,regID,Fx,R1,R2,R3,Z10,Z25,PGA1100)


dzb_if =Cpga(1 );
dzb_sl =Cpga(2 );
nft1   =Cpga(3 );
nft2   =Cpga(4 );
% phi    =Cpga(5 );
% tau    =Cpga(6 );
th1_if =Cpga(7 );
th1_sl =Cpga(8 );
th2_if =Cpga(9 );
th2_sl =Cpga(10);
th3    =Cpga(11);
th4_if =Cpga(12);
th4_sl =Cpga(13);
th5    =Cpga(14);
th6_1  =Cpga(15);
th6_2  =Cpga(16);
th6_3  =Cpga(17);
th6_x1 =Cpga(18);
th6_x2 =Cpga(19);
th6_x3 =Cpga(20);
th6_xc =Cpga(21);
th7    =Cpga(22);
th9_if =Cpga(23);
th9_sl =Cpga(24);
k1     =Cpga(25);
k2     =Cpga(26);
thz1   =Cpga(27);
thz2   =Cpga(28);
thz3   =Cpga(29);
thz4   =Cpga(30);
th11   =Cpga(31);
th12   =Cpga(32);

% Magnitude scaling term
deltaM   = 0.1;
Mref     = 6;
deltaZ   = 1;
th10     = 0;
ZB_if    = 30;
ZB_sl    = 80;

% MAGNITUDE TERM
switch Fs
    case 0, fmag = lh(M,MB,th4_if*(MB - Mref),th4_if,th5,deltaM);
    case 1, fmag = lh(M,MB,th4_sl*(MB - Mref),th4_sl,th5,deltaM);
end

% GEOMETRICAL ATTENUATION TERM
hM = 10.^(nft1+nft2*(M-6));
switch Fs
    case 0, fgeom = (th2_if + th3*M).*log(Rrup+hM);
    case 1, fgeom = (th2_sl + th3*M).*log(Rrup+hM);
end

% DEPTH TERM
Zif_ref = 15;
Zsl_ref = 50;
switch Fs
    case 0, fdepth = lh(Ztor,ZB_if + dzb_if,th9_if*(ZB_if + dzb_if-Zif_ref),th9_if,th10,deltaZ);
    case 1, fdepth = lh(Ztor,ZB_sl + dzb_sl,th9_sl*(ZB_sl + dzb_sl-Zsl_ref),th9_sl,th10,deltaZ);
end

%% SITE TERM
c	= 1.88;
n	= 1.18;
if Vs30<=k1
    fsite = th7*log(Vs30/k1) + k2*(log(PGA1100+c*(Vs30/k1)^n)-log(PGA1100+c));
else
    fsite = (th7+k2*n)*log(Vs30/k1);
end

%% ANAELASTIC ATTENUATION TERM
switch regID
    case 'All'
        fattn = Fx.*(th6_xc+th6_x1*R1+th6_x2*R2+th6_x3*R3)+(1-Fx).*(th6_1*R1+th6_2*R2+th6_3*R3);
    case 'JP'
        fattn = Fx.*(th6_xc+th6_x1*R1+th6_x2*R2+th6_x3*R3)+(1-Fx).*(th6_1*R1+th6_2*R2+th6_3*R3);
    case {'CAM','SA'}
        fattn = Fx.*(th6_xc+th6_x1*R1+th6_x2*R2)+(1-Fx).*(th6_1*R1+th6_2*R2);
    case {'AK','CASC','NZ','TW'}
        fattn = th6_2*Rrup;
    otherwise
        fattn=0;
end

%% BASIN DEPTH TERM
aux    = (log(Vs30)-thz3)/thz4;
lnZref = thz1+(thz2-thz1)*exp(aux)/(1+exp(aux));

switch regID
    case 'CASC'
        delta_lnZ = log(Z25)-lnZref;
    case {'NZ','TW'}
        delta_lnZ = log(Z10)-lnZref;
    case 'JP'
        delta_lnZ = log(Z25)-lnZref;
    otherwise
        th11      = 0;
        th12      = 0;
        delta_lnZ = 0;
end

if Vs30==1100
    fbasin = 0;
else
    fbasin = th11 + th12 *delta_lnZ;
end

%% BASE COEFFICIENT
lnPGA = (1-Fs)*th1_if + Fs*th1_sl + fmag +fgeom + fdepth + fattn + fsite + fbasin;

function y=lh(x,x0,a,b0,b1,delta)

y=a+b0*(x-x0)+(b1-b0)*delta*log(1+exp((x-x0)/delta));

