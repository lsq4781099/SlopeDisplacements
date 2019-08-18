function [lny, sigma,tau,sig] = I_2014_nga(T,M,rrup,mechanism,Vs30)

if T<0 || T> 10
    lny   = nan(size(M));
    sigma = nan(size(M));
    tau   = nan(size(M));
    sig   = nan(size(M));
    %IM    = IM2str(T);
    %h=warndlg(sprintf('GMPE %s not available for %s',mfilename,IM{1}));
    %uiwait(h);
    return
end

To      = max(T,0.001); 
period = [ 0.01	0.02	0.03	0.05	0.075	0.1	0.15	0.2	0.25	0.3	0.4	0.5	0.75	1	1.5	2	3	4	5	7.5	10	0.001];
T_lo    = max(period(period<=To));
T_hi    = min(period(period>=To));
index   = find(abs((period - T_lo)) < 1e-6); % Identify the period

if T_lo==T_hi
    [lny,sigma] = gmpe(index,M,rrup,mechanism,Vs30);
else
    [lny_lo,sigma_lo] = gmpe(index,  M,rrup,mechanism,Vs30);
    [lny_hi,sigma_hi] = gmpe(index+1,M,rrup,mechanism,Vs30);
    x          = log([T_lo;T_hi]);
    Y_sa       = [lny_lo,lny_hi]';
    Y_sigma    = [sigma_lo,sigma_hi]';
    lny        = interp1(x,Y_sa,log(To))';
    sigma      = interp1(x,Y_sigma,log(To))';
end
tau = sigma*0;
sig = sigma;


function[lny,sigma]=gmpe(index,M,rrup,mechanism,Vs30)

period = [ 0.01	0.02	0.03	0.05	0.075	0.1	0.15	0.2	0.25	0.3	0.4	0.5	0.75	1	1.5	2	3	4	5	7.5	10	0.001];
T = period(index);
if Vs30 > 1200
    Vs30 = 1200;
end

M_SE = max(min(M, 7.5),5);
if T <= 0.05
    sigma = 1.18 + 0.035 * log(0.05) - 0.06 * M_SE;
elseif T >= 3
    sigma = 1.18 + 0.035 * log(3) - 0.06 * M_SE;
else
    sigma = 1.18 + 0.035 * log(T) - 0.06 * M_SE;
end

switch mechanism
    case 'strike-slip',F=0;
    case 'reverse',    F=1;
end

Coef1 = [7.0887 7.1157 7.2087 6.2638 5.9051 7.5791 8.0190 9.2812 9.5804 9.8912 9.5342 9.2142 8.3517 7.0453 5.1307 3.3610 0.1784 -2.4301 -4.3570 -7.8275 -9.2857 7.0887
    0.2058 0.2058 0.2058 0.0625 0.1128 0.0848 0.1713 0.1041 0.0875 0.0003 0.0027 0.0399 0.0689 0.1600 0.2429 0.3966 0.7560 0.9283 1.1209 1.4016 1.5574 0.2058
    2.9935 2.9935 2.9935 2.8664 2.9406 3.0190 2.7871 2.8611 2.8289 2.8423 2.8300 2.8560 2.7544 2.7339 2.6800 2.6837 2.6907 2.5782 2.5468 2.4478 2.3922 2.9935
    -0.2287 -0.2287 -0.2287 -0.2418 -0.2513 -0.2516 -0.2236 -0.2229 -0.2200 -0.2284 -0.2318 -0.2337 -0.2392 -0.2398 -0.2417 -0.2450 -0.2389 -0.2514 -0.2541 -0.2593 -0.2586 -0.2287];
Coef2 = [ 9.0138 9.0408 9.1338 7.9837 7.7560 9.4252 9.6242 11.1300 11.3629 11.7818 11.6097 11.4484 10.9065 9.8565 8.3363 6.8656 4.1178 1.8102 0.0977 -3.0563 -4.4387 9.0138
    -0.0794 -0.0794 -0.0794 -0.1923 -0.1614 -0.1887 -0.0665 -0.1698 -0.1766 -0.2798 -0.3048 -0.2911 -0.3097 -0.2565 -0.2320 -0.1226 0.1724 0.3001 0.4609 0.6948 0.8393 -0.0794
    2.9935 2.9935 2.9935 2.7995 2.8143 2.8131 2.4091 2.4938 2.3773 2.3772 2.3413 2.3477 2.2042 2.1493 2.0408 2.0013 1.9408 1.7763 1.7030 1.5212 1.4195 2.9935
    -0.2287 -0.2287 -0.2287 -0.2319 -0.2326 -0.2211 -0.1676 -0.1685 -0.1531 -0.1595 -0.1594 -0.1584 -0.1577 -0.1532 -0.1470 -0.1439 -0.1278 -0.1326 -0.1291 -0.1220 -0.1145 -0.2287];

a3	  = [ 0.0589 0.0589 0.0589 0.0417 0.0527 0.0442 0.0329 0.0188 0.0095 -0.0039 -0.0133 -0.0224 -0.0267 -0.0198 -0.0367 -0.0291 -0.0214 -0.0240 -0.0202 -0.0219 -0.0035 0.0589];
zeta  = [ -0.854 -0.854 -0.854 -0.631 -0.591 -0.757 -0.911 -0.998 -1.042 -1.030 -1.019 -1.023 -1.056 -1.009 -0.898 -0.851 -0.761 -0.675 -0.629 -0.531 -0.586 -0.854];
gamma = [ -0.0027 -0.0027 -0.0027 -0.0061 -0.0056 -0.0042 -0.0046 -0.0030 -0.0028 -0.0029 -0.0028 -0.0021 -0.0029 -0.0032 -0.0033 -0.0032 -0.0031 -0.0051 -0.0059 -0.0057 -0.0061 -0.0027];
phi   = [ 0.08 0.08 0.08 0.08 0.08 0.08 0.08 0.08 0.08 0.08 0.08 0.08 0.08 0.06 0.04 0.02 0.02 0 0 0 0 0.08];

C1    = Coef1(:,index);
C2    = Coef2(:,index);
a1    = (M<=6.75)*C1(1)+and(6.75<M,M<=8.5)*C2(1);
a2    = (M<=6.75)*C1(2)+and(6.75<M,M<=8.5)*C2(2);
b1    = (M<=6.75)*C1(3)+and(6.75<M,M<=8.5)*C2(3);
b2    = (M<=6.75)*C1(4)+and(6.75<M,M<=8.5)*C2(4);
a3    = a3(index);
zeta  = zeta(index);
gamma = gamma(index);
phi   = phi(index);
lny   = a1 + a2 .* M + a3 .* (8.5 - M).^2 - (b1 + b2 .* M).* log(rrup + 10) + zeta .* log(Vs30) + gamma .* rrup + phi*F;

