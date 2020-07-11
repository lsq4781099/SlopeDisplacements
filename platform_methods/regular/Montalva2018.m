function [lny,sigma,tau,phi]=Montalva2018(To,M,Rrup,Rhyp,Zhyp,VS30,f0,mechanism)

st = dbstack;
[isadmisible,units] = isIMadmisible(To,st(1).name,[0 10],[nan nan],[nan nan],[nan nan]);
if isadmisible==0
    lny   = nan(size(M));
    sigma = nan(size(M));
    tau   = nan(size(M));
    phi   = nan(size(M));
    return
end
To      = max(To,0.01); %PGA is associated to To=0.01;
period  = [0.01,0.02,0.05,0.075,0.1,0.2,0.25,0.4,0.5,0.6,0.75,1,1.5,2,2.5,3,4,5,6,7.5,10];
T_lo    = max(period(period<=To));
T_hi    = min(period(period>=To));
index   = find(abs((period - T_lo)) < 1e-6); % Identify the period

switch mechanism
    case 'interface', R=Rrup; fevent=0;
    case 'intraslab', R=Rhyp; fevent=1;
end

if f0==-999
    fpeak=0;
else
    fpeak=1;
end

Arock=exp(gmpe(1,M,R,Zhyp,1000,f0,fevent,fpeak));

if T_lo==T_hi
    [lny,sigma,tau] = gmpe(index,M,R,Zhyp,VS30,f0,fevent,fpeak,Arock);
    phi = sqrt(sigma.^2-tau.^2);
else
    [lny_lo,sigma_lo,tau_lo] = gmpe(index,  M,R,Zhyp,VS30,f0,fevent,fpeak,Arock);
    [lny_hi,sigma_hi,tau_hi] = gmpe(index+1,M,R,Zhyp,VS30,f0,fevent,fpeak,Arock);
    x          = log([T_lo;T_hi]);
    Y_sa       = [lny_lo,lny_hi]';
    Y_sigma    = [sigma_lo,sigma_hi]';
    Y_tau      = [tau_lo,tau_hi]';
    lny        = interp1(x,Y_sa,log(To))';
    sigma      = interp1(x,Y_sigma,log(To))';
    tau        = interp1(x,Y_tau,log(To))';
    phi        = sqrt(sigma.^2-tau.^2);
end

% unit convertion
lny  = lny+log(units);

function[lny,sigma,tau]=gmpe(index,M,R,Zhyp,VS30,f0,fevent,fpeak,Arock)

%Coeficientes
Coef = [9.5796	1.694	-2.0552	-0.6719	0.0481	-0.0029	0.0086	4.3234	-0.1836	1.5342	-0.2227	0.0681	0.5496	0.3762	0.3205	0.7392
7.2827	0.4879	-1.6088	-0.4773	0.2456	-0.0042	0.0114	3.1605	-0.181	-0.2027	-0.2572	0.063	0.5679	0.387	0.3189	0.7576
2.3273	2.2898	-0.5017	-1.0245	-0.1336	-0.0115	0.0086	6.2085	-0.1462	0.784	0.0553	0.0301	0.6166	0.3355	0.3289	0.7752
2.7371	3.4236	-0.5182	-1.2203	-0.3276	-0.0118	0.0061	7.4944	0.0153	2.0237	-0.0005	-0.0104	0.6663	0.3511	0.3328	0.8234
12.3513	6.2946	-2.3478	-1.3019	-0.7962	-0.0067	0.0033	8.1111	0.0092	6.89	0.0542	0.0095	0.6543	0.387	0.3318	0.8294
12.0135	3.0978	-2.3538	-0.6445	-0.1895	-0.0031	0.0039	4.5146	0.0012	3.4707	-0.4038	0.0493	0.6146	0.4375	0.3353	0.8255
5.8869	2.2804	-1.1626	-0.7092	-0.0846	-0.0065	0.0004	4.6583	-0.1023	1.5693	-0.4341	0.042	0.5437	0.4618	0.3381	0.7894
7.7585	2.6098	-1.6313	-0.6409	-0.1235	-0.0029	0.0008	3.8693	-0.3932	2.3911	-0.0972	0.023	0.5208	0.5546	0.356	0.8399
11.0256	3.4167	-2.3273	-0.673	-0.2398	0.0001	-0.0006	3.915	-0.6166	3.964	0.0253	0.0321	0.4659	0.4692	0.3629	0.7542
14.3245	4.0743	-3.0383	-0.8353	-0.3266	0.003	0.0018	4.56	-0.8733	5.2928	0.3383	0.0236	0.405	0.3975	0.3765	0.681
13.0566	3.9789	-2.8202	-0.8006	-0.2936	0.0018	0.0017	4.2968	-1.0607	4.7228	0.5661	0.0133	0.4103	0.3553	0.3646	0.6538
10.1588	3.5206	-2.3521	-0.8777	-0.2049	0.0015	0.0003	4.6197	-1.0904	3.7595	0.673	-0.0363	0.5375	0.3942	0.365	0.76
6.575	4.2188	-1.648	-1.1342	-0.3511	-0.0013	0.0018	5.9613	-1.0085	3.9026	0.666	-0.0634	0.5609	0.4161	0.3592	0.7853
12.4348	7.6973	-2.8442	-1.5506	-0.9359	0.0017	0.001	8.1448	-0.5326	8.7793	0.4262	-0.0917	0.5679	0.3963	0.3474	0.7748
16.4865	7.1936	-3.7101	-1.4618	-0.7409	0.0072	-0.0037	7.8425	-0.2571	8.3016	0.2749	-0.0806	0.5607	0.3782	0.3526	0.7628
19.2618	7.6039	-4.3416	-1.1854	-0.7872	0.0098	-0.0035	6.4895	-0.0947	9.4182	0.0654	-0.0708	0.5248	0.3995	0.3446	0.7442
14.9964	6.1593	-3.6399	-0.5839	-0.5437	0.007	-0.0035	3.4696	-0.1265	7.3303	0.25	-0.0515	0.4933	0.3141	0.3246	0.6689
18.7614	7.4247	-4.4797	-0.6542	-0.7506	0.0101	-0.0019	3.8682	-0.2322	9.5312	0.3209	-0.0268	0.4927	0.3454	0.2872	0.6667
20.7315	7.9061	-4.9686	-0.9894	-0.8281	0.0117	-0.0005	5.7567	-0.3653	10.8509	0.2682	0.0409	0.4767	0.2928	0.277	0.6243
18.6702	7.832	-4.6018	-1.0741	-0.8287	0.0097	-0.0036	6.2898	-0.3176	10.6605	0.1471	0.1168	0.4563	0.3596	0.2881	0.6484
16.4001	8.1968	-4.3047	-1.2249	-0.9383	0.0089	-0.005	7.1431	-0.323	11.5219	0.3459	0.0484	0.4987	0.249	0.2785	0.6231];


a1     = Coef(index,1);
a2     = Coef(index,2);
a4     = Coef(index,3);
a5     = Coef(index,4);
a6     = Coef(index,5);
a8     = Coef(index,6);
a9     = Coef(index,7);
a10    = Coef(index,8);
a11    = Coef(index,9);
a12    = Coef(index,10);
a13    = Coef(index,11);
a14    = Coef(index,12);
tau    = Coef(index,13);
sigma  = Coef(index,16);

C1   = 7.4;
C2   = 20;
C3   = 0.95;
Vlin = 725;
flin = 5.5;
c    = 4.0;
n1   = 1.8;
n2   = 3.0;


% source term
ind         = M>=C1;
Ffuente     = a2*(M-C1);
Ffuente(ind)= a12*(M(ind)-C1);

% Path Term
Ftrayectoria = (a4+a5*fevent+a6*(M-C1)).*log(R+C2*exp(C3*(M-5)))+a8*R;

% Event term
Zh      = min(Zhyp,120);
Fevento = (a10+a9*(Zh-60))*fevent;

% Site term
VS  = min(1000, VS30); % (equation (8))
f0S = min(6.0,f0);
if VS30 < Vlin
    Fsitio = a11*log(VS/Vlin)+a14*fpeak*log(f0S/flin)-a13*log(Arock+c)+a13*(log(Arock+c*((VS/Vlin)^n1+fpeak*(f0S/flin)^n2)));
else
    Fsitio = a11*log(VS/Vlin)+a14*fpeak*log(f0S/flin)+a13*(n1*log(VS/Vlin)+n2*fpeak*log(f0S/flin));
end


lny = a1 + Ffuente + Ftrayectoria + Fevento + Fsitio;
