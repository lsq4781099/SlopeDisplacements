function [IM,Rmetric,Residuals]=mGMPE_info(model)
%
% [IM,Rmetric,Residuals]=mGMPE_info(model)
%
% Function that retuns information about GMPEs
%
% IM = Intensity Measure Id.
%    =  0  for PGA      (Peak Ground Acceleration)
%    = -1  for PGV      (Peak Ground Velocity)
%    = -2  for PGD      (Permanent Ground Displacement)
%    = -3  for Duration (Duration)
%    = -4  for CAV      (Cummulative Absolute Velocity)
%    = -5  for AI       (Arias Intensity)
%    = -6  for VGI      (Incremental Ground Velocity)
%    = To    for Sa     (Pseudoacceleration    at period To , To>0)
%    = To+i  for Sv     (Pseudo Celocity       at period To , To>0)
%    = To+2i for Sd     (Spectral Displacement at period To , To>0)
%    = To+3i for H/V    (H/V ratio             at period To , To>0)

% Rmetric is a vector of the distance metrics used by the GMPE (1 or 0)
% if Rmetric(i)=1 means that the i-th metric will be computed, and 0 will
% not be computed. The metris available are:
%
% Rmetric    = [Rrup Rhyp Rjb Repi Zseis Rx Ry0 zhyp ztor zbor zbot]
% 1.-  rrup  = closest diatance from site to rupture plane
% 2.-  rhyp  = distance from the site to the hypocenter
% 3.-  rjb   = Joyner-Boore distance, i.e., closest distance from site to
%              surface projection of rupture area
% 4.-  repi  = epicentral distance
% 5.-  rseis = Shortest distance between the recording site and the presumed
%              zone of seismogenic rupture on the fault
% 6.-  rx    = Horizontal distance from top of rupture measured perpendicular
%              to fault strike
% 7.-  ry0   = Horizontal distance off the end of the rupture measured
%              parallel to strike
% 8.-  zhyp  = depth of hypocenter
% 9.-  ztor  = Depth to top of coseismic rupture (km)
% 10.- zbor  = Depth to the bottom of the rupture plane (km)
% 11.- zbot  = The depth to the bottom of the seismogenic crust (km)

%
% Residuals = Probability distribution of IM rsiduals, e.g., lognormal
switch lower(model)
    case 'youngs1997',                 IM = [0.01 0.075 0.1 0.2 0.3 0.4 0.5 0.75 1 1.5 2 3];
    case 'atkinsonboore2003',          IM = [0.01 0.04 0.1 0.2 0.4 1 2 1/0.33];
    case 'zhao2006',                   IM = [0.02 0.05 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 1.25 1.5 2 2.5 3 4 5];
    case 'mcverry2006',                IM = [0.001 0.075 0.1 0.2 0.3 0.4 0.5 0.75 1 1.5 2 3];
    case 'bchydro2012',                IM = [0.01 0.02 0.05 0.075 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.6 0.75 1 1.5 2 2.5 3 4 5 6 7.5 10];
    case 'bchydro2018',                IM = [0.01 0.02 0.03 0.05 0.075 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.6 0.75 1 1.5 2 2.5 3 4 5 6 7.5 10];
    case 'kuehn2020',                  IM = [-1 0 0.01 0.02 0.03 0.05 0.075 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.75 1 1.5 2 3 4 5 7.5 10];
    case 'parker2020',                 IM = [-1 0.001 0.01 0.02 0.03 0.05 0.075 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.75 1 1.5 2 3 4 5 7.5 10];
    case 'arteta2018',                 IM = [0.01 0.02 0.05 0.075 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.6 0.75 1 1.5 2 2.5 3 4 5 6 7.5 10];
    case 'montalvabastias2017',        IM = [0.01 0.02 0.05 0.075 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.6 0.75 1 1.5 2 2.5 3 4 5 6 7.5 10];
    case 'montalvabastias2017hq',      IM = [0.01 0.02 0.05 0.075 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.6 0.75 1 1.5 2 2.5 3 4 5 6 7.5 10];
    case 'montalva2018',               IM = [0.01 0.02 0.05 0.075 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.6 0.75 1 1.5 2 2.5 3 4 5 6 7.5 10];
    case 'siberrisk2019',              IM = [-1 0.01 0.02 0.05 0.075 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.6 0.75 1 1.5 2 2.5 3 4 5 6 7.5 10];    
    case 'idini2016',                  IM = [0.001 0.01 0.02 0.03 0.05 0.07 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.75 1 1.5 2 3 4 5 7.5 10];
    case 'contrerasboroschek2012',     IM = [0.01 0.04 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 1.1 1.2 1.3 1.4 1.5 2];
    case 'garcia2005',                 IM = [-1 0 1/25 1/20 1/13.33 1/10 1/5 1/3.33 1/2.5 1/2 1/1.33 1/1 1/0.67 1/0.5 1/0.33 1/0.25 1/0.2];
    case 'jaimes2006',                 IM = [0.01 0.2 0.4 0.6 0.8 1 1.2 1.4 1.6 1.8 2 2.2 2.4 2.6 2.8 3 3.2 3.4 3.6 3.8 4 4.2 4.4 4.6 4.8 5 5.2 5.4 5.6 5.8 6];
    case 'jaimes2015',                 IM = [-1 0.01 0.2 0.4 0.6 0.8 1 1.2 1.4 1.6 1.8 2 2.2 2.4 2.6 2.8 3 3.2 3.4 3.6 3.8 4 4.2 4.4 4.6 4.8 5];
    case 'jaimes2016',                 IM = [-1 0.001 0.01 0.02 0.03 0.04 0.05 0.08 0.1 0.12 0.15 0.17 0.2 0.25 0.3 0.4 0.5 0.75 1 1.5 2 3 4 5 7.5 10];
    case 'garciajaimes2017',           IM = [-1 0.001 0.01 0.02 0.04 0.06 0.08 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2 2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8 2.9 3 4 5];
    case 'garciajaimes2017hv',         IM =   3i+[0.01 0.02 0.04 0.06 0.08 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2 2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8 2.9 3 4 5];
    case 'bernal2014',                 IM = [0.001 0.05 0.1 0.15 0.3 0.5 1 1.5 2 2.5 3 7 8];
    case 'sadigh1997',                 IM = [0.02 0.075 0.10 0.20 0.30 0.40 0.50 0.75 1 1.50 2 3 4];
    case 'i2008',                      IM = [0.01 0.02 0.03 0.05 0.075 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.75 1 1.5 2 3 4 5 7.5 10];
    case 'cy2008',                     IM = [-1 0.001 0.01 0.02 0.03 0.04 0.05 0.075 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.75 1 1.5 2 3 4 5 7.5 10];
    case 'ba2008',                     IM = [-1 0.001 0.01 0.02 0.03 0.05 0.075 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.75 1 1.5 2 3 4 5 7.5 10];
    case 'cb2008',                     IM = [-2 -1 0.001 0.01 0.02 0.03 0.05 0.075 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.75 1 1.5 2 3 4 5 7.5 10];
    case 'as2008',                     IM = [-1 0 0.01 0.02 0.03 0.04 0.05 0.075 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.75 1 1.5 2 3 4 5 7.5 10];
    case 'as1997h',                    IM = [0.01 0.02 0.03 0.04 0.05 0.06 0.075 0.09 0.1 0.12 0.15 0.17 0.2 0.24 0.3 0.36 0.4 0.46 0.5 0.6 0.75 0.85 1 1.5 2 3 4 5];
    case 'i2014',                      IM = [0.001 0.01 0.02 0.03 0.05 0.075 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.75 1 1.5 2 3 4 5 7.5 10];
    case 'cy2014',                     IM = [-1 0.001 0.01 0.02 0.03 0.04 0.05 0.075 0.1 0.12 0.15 0.17 0.2 0.25 0.3 0.4 0.5 0.75 1 1.5 2 3 4 5 7.5 10];
    case 'bssa2014',                   IM = [-1 0.001 0.01 0.02 0.03 0.05 0.075 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.75 1 1.5 2 3 4 5 7.5 10];
    case 'cb2014',                     IM = [-1 0.001 0.010 0.020 0.030 0.050 0.075 0.10 0.15 0.20 0.25 0.30 0.40 0.50 0.75 1 1.5 2 3 4 5 7.5 10];
    case 'ask2014',                    IM = [-1 0.001 0.01 0.02 0.03 0.05 0.075 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.75 1 1.5 2 3 4 5 6 7.5 10];
    case 'akkarboomer2007',            IM = [0 [0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1 1.05 1.1 1.15 1.2 1.25 1.3 1.35 1.4 1.45 1.5 1.55 1.6 1.65 1.7 1.75 1.8 1.85 1.9 1.95 2 2.05 2.1 2.15 2.2 2.25 2.3 2.35 2.4 2.45 2.5 2.55 2.6 2.65 2.7 2.75 2.8 2.85 2.9 2.95 3 3.05 3.1 3.15 3.2 3.25 3.3 3.35 3.4 3.45 3.5 3.55 3.6 3.65 3.7 3.75 3.8 3.85 3.9 3.95 4]+2*1i];
    case 'akkarboomer2010',            IM = [-1 0:0.05:3];
    case 'arroyo2010',                 IM = [0 0.04 0.045 0.05 0.055 0.06 0.065 0.07 0.075 0.08 0.085 0.09 0.095 0.1 0.12 0.14 0.16 0.18 0.2 0.22 0.24 0.26 0.28 0.3 0.32 0.34 0.36 0.38 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2 2.5 3 3.5 4 4.5 5];
    case 'bindi2011',                  IM = [-1 0 0.04 0.07 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.6 0.7 0.8 0.9 1 1.25 1.5 1.75 2];        
    case 'kanno2006',                  IM = [-1 0 0.05 0.06 0.07 0.08 0.09 0.1 0.11 0.12 0.13 0.15 0.17 0.2 0.22 0.25 0.3 0.35 0.4 0.45 0.5 0.6 0.7 0.8 0.9 1 1.1 1.2 1.3 1.5 1.7 2 2.2 2.5 3 3.5 4 4.5 5];
    case 'cauzzi2015',                 IM = [-1 round([0.01:0.01:0.09 0.1:0.05:10]*100)/100 round([0.01:0.01:0.09 0.1:0.05:10]*100)/100+2i];
    case 'dw12',                       IM = -4;
    case 'fg15',                       IM = [-5 -4];
    case 'tba03',                      IM = -5;
    case 'bu17',                       IM = [-6 -5 -4];
    case 'cb10',                       IM = -4;
    case 'cb11',                       IM = -4;
    case 'cb19',                       IM = [-5 -4];
    case 'km06',                       IM = -4;
    case 'medianpce_bchydro',          IM = [0.01 0.02 0.05 0.075 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.6 0.75 1 1.5 2 2.5 3];
    case 'medianpce_nga',              IM = [0.001 0.02 0.03 0.05 0.075 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.75 1 1.5 2 3 4 5 6 7 8 10];

    case 'udm',                        IM = [];
    case 'macedo2019',                 IM = -5;
    case 'macedo2020',                 IM = -4;
    case 'pce_nga',                    IM = [0.001 0.02 0.03 0.05 0.075 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.75 1 1.5 2 3 4 5 6 7 8 10];
    case 'pce_bchydro',                IM = [0.01 0.02 0.05 0.075 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.6 0.75 1 1.5 2 2.5 3];
    case 'franky',                     IM = [];
   
end

if nargout>=2
    switch lower(model)
        case 'youngs1997',                 Rmetric=[1 0 0 0 0 0 0 1 0 0 0];  Residuals = 'lognormal';
        case 'atkinsonboore2003',          Rmetric=[1 0 0 0 0 0 0 1 0 0 0];  Residuals = 'lognormal';
        case 'zhao2006',                   Rmetric=[1 0 0 0 0 0 0 1 0 0 0];  Residuals = 'lognormal';
        case 'mcverry2006',                Rmetric=[1 0 0 0 0 0 0 1 0 0 0];  Residuals = 'lognormal';
        case 'contrerasboroschek2012',     Rmetric=[1 0 0 0 0 0 0 1 0 0 0];  Residuals = 'lognormal';
        case 'bchydro2012',                Rmetric=[1 1 0 0 0 0 0 1 0 0 0];  Residuals = 'lognormal';
        case 'bchydro2018',                Rmetric=[1 0 0 0 0 0 0 0 1 0 0];  Residuals = 'lognormal';
        case 'kuehn2020',                  Rmetric=[1 0 0 0 0 0 0 0 1 0 0];  Residuals = 'lognormal';
        case 'parker2020',                 Rmetric=[1 0 0 0 0 0 0 1 0 0 0];  Residuals = 'lognormal';
        case 'arteta2018',                 Rmetric=[0 1 0 0 0 0 0 0 0 0 0];  Residuals = 'lognormal';
        case 'idini2016',                  Rmetric=[1 1 0 0 0 0 0 1 0 0 0];  Residuals = 'lognormal';
        case 'montalvabastias2017',        Rmetric=[1 1 0 0 0 0 0 1 0 0 0];  Residuals = 'lognormal';
        case 'montalvabastias2017hq',      Rmetric=[1 1 0 0 0 0 0 1 0 0 0];  Residuals = 'lognormal';
        case 'montalva2018',               Rmetric=[1 1 0 0 0 0 0 1 0 0 0];  Residuals = 'lognormal';
        case 'siberrisk2019',              Rmetric=[1 1 0 0 0 0 0 1 0 0 0];  Residuals = 'lognormal';
        case 'garcia2005',                 Rmetric=[1 1 0 0 0 0 0 1 0 0 0];  Residuals = 'lognormal';
        case 'jaimes2006',                 Rmetric=[1 0 0 0 0 0 0 1 0 0 0];  Residuals = 'lognormal';
        case 'jaimes2015',                 Rmetric=[1 0 0 0 0 0 0 1 0 0 0];  Residuals = 'lognormal';
        case 'jaimes2016',                 Rmetric=[1 0 0 0 0 0 0 1 0 0 0];  Residuals = 'lognormal';
        case 'garciajaimes2017',           Rmetric=[1 0 0 0 0 0 0 1 0 0 0];  Residuals = 'lognormal';
        case 'garciajaimes2017hv',         Rmetric=[1 0 0 0 0 0 0 1 0 0 0];  Residuals = 'lognormal';
        case 'bernal2014',                 Rmetric=[1 0 0 0 0 0 0 1 0 0 0];  Residuals = 'lognormal';
        case 'sadigh1997',                 Rmetric=[1 0 0 0 0 0 0 1 0 0 0];  Residuals = 'lognormal';
        case 'i2008',                      Rmetric=[1 0 0 0 0 0 0 1 0 0 0];  Residuals = 'lognormal';
        case 'cy2008',                     Rmetric=[1 0 1 0 0 1 0 0 1 0 0];  Residuals = 'lognormal';
        case 'ba2008',                     Rmetric=[0 0 1 0 0 0 0 0 0 0 0];  Residuals = 'lognormal';
        case 'cb2008',                     Rmetric=[1 0 1 0 0 0 0 0 1 0 0];  Residuals = 'lognormal';
        case 'as2008',                     Rmetric=[1 0 1 0 0 1 0 0 1 0 0];  Residuals = 'lognormal';
        case 'as1997h',                    Rmetric=[1 0 0 0 0 0 0 1 0 0 0];  Residuals = 'lognormal';
        case 'i2014',                      Rmetric=[1 0 0 0 0 0 0 0 0 0 0];  Residuals = 'lognormal';
        case 'cy2014',                     Rmetric=[1 0 1 0 0 1 0 0 1 0 0];  Residuals = 'lognormal';
        case 'bssa2014',                   Rmetric=[0 0 1 0 0 0 0 0 0 0 0];  Residuals = 'lognormal';
        case 'cb2014',                     Rmetric=[1 0 1 0 0 1 0 1 1 0 1];  Residuals = 'lognormal';
        case 'ask2014',                    Rmetric=[1 0 1 0 0 1 1 0 1 0 0];  Residuals = 'lognormal';
        case 'akkarboomer2007',            Rmetric=[0 0 1 0 0 0 0 0 0 0 0];  Residuals = 'lognormal';
        case 'akkarboomer2010',            Rmetric=[0 0 1 0 0 0 0 0 0 0 0];  Residuals = 'lognormal';
        case 'arroyo2010',                 Rmetric=[1 0 0 0 0 0 0 0 0 0 0];  Residuals = 'lognormal';
        case 'bindi2011',                  Rmetric=[0 0 1 0 0 0 0 0 0 0 0];  Residuals = 'lognormal';
        case 'kanno2006',                  Rmetric=[1 0 0 0 0 0 0 1 0 0 0];  Residuals = 'lognormal';
        case 'cauzzi2015',                 Rmetric=[1 1 0 0 0 0 0 0 0 0 0];  Residuals = 'lognormal';
        case 'dw12',                       Rmetric=[1 0 0 0 0 0 0 0 0 0 0];  Residuals = 'lognormal';
        case 'fg15',                       Rmetric=[1 0 0 0 0 0 0 1 0 0 0];  Residuals = 'lognormal';
        case 'tba03',                      Rmetric=[1 0 0 0 0 0 0 0 0 0 0];  Residuals = 'lognormal';
        case 'bu17',                       Rmetric=[1 0 0 0 0 0 0 1 0 0 0];  Residuals = 'lognormal';
        case 'cb10',                       Rmetric=[1 0 1 0 0 0 0 0 1 0 0];  Residuals = 'lognormal';
        case 'cb11',                       Rmetric=[1 0 1 0 0 0 0 0 1 0 0];  Residuals = 'lognormal';
        case 'cb19',                       Rmetric=[1 0 1 0 0 1 0 1 1 0 1];  Residuals = 'lognormal';
        case 'km06',                       Rmetric=[1 0 0 0 0 0 0 0 0 0 0];  Residuals = 'lognormal';
        case 'medianpce_bchydro',          Rmetric=[1 0 0 0 0 0 0 0 0 0 0];  Residuals = 'lognormal';
        case 'medianpce_nga',              Rmetric=[1 0 0 0 0 0 0 0 0 0 0];  Residuals = 'lognormal';
        case 'udm',                        Rmetric=[];                       Residuals = '';
        case 'macedo2019',                 Rmetric=[1 0 0 0 0 0 0 0 0 0 0];  Residuals = 'lognormal';    
        case 'macedo2020',                 Rmetric=[1 0 0 0 0 0 0 0 0 0 0];  Residuals = 'lognormal';
        case 'pce_nga',                    Rmetric=[1 0 0 0 0 0 0 0 0 0 0];  Residuals = 'lognormal';
        case 'pce_bchydro',                Rmetric=[1 0 0 0 0 0 0 0 0 0 0];  Residuals = 'lognormal';    
        case 'franky',                     Rmetric=[];                       Residuals = '';            
    end
end

