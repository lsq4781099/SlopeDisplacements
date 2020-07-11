function ME = pshatoolbox_methods(group,varargin)

% regular = typical GMPEs
% cond    = conditional GMPE
% udm     = user defined models, GMPEs defined in user specified m-files
% pce     = polynomial chaos expansion GMPE

% ground motion models (gmpe)
if group ==1
    i = 0;
    ME(1:49,1)=struct('label',[],'str',[],'func',[],'type',[],'mech',[],'ref',[]);
    % regular GMMs
    i=i+1; ME(i).label = 'Youngs et al. 1997';           ME(i).func = @Youngs1997;                 ME(i).mech=[1 1 0]; ME(i).type='regular';  ME(i).ref='https://doi.org/10.1785/gssrl.68.1.58';
    i=i+1; ME(i).label = 'Atkinson & Boore, 2003';       ME(i).func = @AtkinsonBoore2003;          ME(i).mech=[1 1 0]; ME(i).type='regular';  ME(i).ref='https://doi.org/10.1785/0120020156';
    i=i+1; ME(i).label = 'Zhao et al. 2006';             ME(i).func = @Zhao2006;                   ME(i).mech=[1 1 1]; ME(i).type='regular';  ME(i).ref='https://doi.org/10.1785/0120050122';
    i=i+1; ME(i).label = 'McVerry et al. 2006';          ME(i).func = @Mcverry2006;                ME(i).mech=[1 1 1]; ME(i).type='regular';  ME(i).ref='http://www.nzsee.org.nz/db/Bulletin/Archive/39(1)0001.pdf';
    i=i+1; ME(i).label = 'Boroschek et al. 2012';        ME(i).func = @ContrerasBoroschek2012;     ME(i).mech=[1 0 0]; ME(i).type='regular';  ME(i).ref='https://nisee.berkeley.edu/elibrary/files/documents/elib/www/documents/201204/PISELL/boroschek-maule-eq.pdf';
    i=i+1; ME(i).label = 'Abrahamson et al. 2016';       ME(i).func = @BCHydro2012;                ME(i).mech=[1 1 0]; ME(i).type='regular';  ME(i).ref='https://doi.org/10.1193/051712EQS188MR';
    i=i+1; ME(i).label = 'Abrahamson et al. 2018';       ME(i).func = @BCHydro2018;                ME(i).mech=[1 1 0]; ME(i).type='regular';  ME(i).ref='https://peer.berkeley.edu/sites/default/files/2018_02_abrahamson_9.10.18.pdf';
    i=i+1; ME(i).label = 'Kuehn et al. 2020';            ME(i).func = @Kuehn2020;                  ME(i).mech=[1 1 0]; ME(i).type='regular';  ME(i).ref='https://peer.berkeley.edu/sites/default/files/2020_02_final_3.17.2020.pdf';
    i=i+1; ME(i).label = 'Parker et al. 2020';           ME(i).func = @Parker2020;                 ME(i).mech=[1 1 0]; ME(i).type='regular';  ME(i).ref='https://peer.berkeley.edu/sites/default/files/2020_02_final_3.17.2020.pdf';
    i=i+1; ME(i).label = 'Arteta et al. 2018';           ME(i).func = @Arteta2018;                 ME(i).mech=[0 1 0]; ME(i).type='regular';  ME(i).ref='https://doi.org/10.1193/102116EQS176M';
    i=i+1; ME(i).label = 'Idini et al. 2016';            ME(i).func = @Idini2016;                  ME(i).mech=[1 1 0]; ME(i).type='regular';  ME(i).ref='https://doi.org/10.1007/s10518-016-0050-1';
    i=i+1; ME(i).label = 'Montalva et al. 2017';         ME(i).func = @MontalvaBastias2017;        ME(i).mech=[1 1 0]; ME(i).type='regular';  ME(i).ref='https://doi.org/10.1785/0120160221';
    i=i+1; ME(i).label = 'Montalva et al. 2017 (HQ) ';   ME(i).func = @MontalvaBastias2017HQ;      ME(i).mech=[1 1 0]; ME(i).type='regular';  ME(i).ref='https://doi.org/10.1785/0120160221';
    i=i+1; ME(i).label = 'Montalva et al. 2018';         ME(i).func = @Montalva2018;               ME(i).mech=[1 1 0]; ME(i).type='regular';  ME(i).ref='';
    i=i+1; ME(i).label = 'SIBER-RISK 2019';              ME(i).func = @SiberRisk2019;              ME(i).mech=[1 1 0]; ME(i).type='regular';  ME(i).ref='https://doi.org/10.1177/8755293019891723';
    i=i+1; ME(i).label = 'Garcia et al. 2005';           ME(i).func = @Garcia2005;                 ME(i).mech=[0 1 0]; ME(i).type='regular';  ME(i).ref='https://doi.org/10.1785/0120050072';
    i=i+1; ME(i).label = 'Jaimes et al. 2006';           ME(i).func = @Jaimes2006;                 ME(i).mech=[1 0 0]; ME(i).type='regular';  ME(i).ref='https://doi.org/10.1080/13632460609350622';
    i=i+1; ME(i).label = 'Jaimes et al. 2015';           ME(i).func = @Jaimes2015;                 ME(i).mech=[0 1 0]; ME(i).type='regular';  ME(i).ref='https://doi.org/10.1080/13632469.2015.1025926';
    i=i+1; ME(i).label = 'Jaimes et al. 2016';           ME(i).func = @Jaimes2016;                 ME(i).mech=[0 0 1]; ME(i).type='regular';  ME(i).ref='https://doi.org/10.1785/0120150283';
    i=i+1; ME(i).label = 'Garcia-Soto Jaimes 2017';      ME(i).func = @GarciaJaimes2017;           ME(i).mech=[1 0 0]; ME(i).type='regular';  ME(i).ref='https://doi.org/10.1785/0120160273';
    i=i+1; ME(i).label = 'Garcia-Soto Jaimes 2017 (H/V)';ME(i).func = @GarciaJaimes2017HV;         ME(i).mech=[1 0 0]; ME(i).type='regular';  ME(i).ref='https://doi.org/10.1785/0120160273';
    i=i+1; ME(i).label = 'Bernal et al. 2014';           ME(i).func = @Bernal2014;                 ME(i).mech=[1 1 0]; ME(i).type='regular';  ME(i).ref='https://doi.org/10.13140/2.1.2693.6641';
    i=i+1; ME(i).label = 'Sadigh et al. 1997';           ME(i).func = @Sadigh1997;                 ME(i).mech=[0 0 1]; ME(i).type='regular';  ME(i).ref='https://doi.org/10.1785/gssrl.68.1.180';
    i=i+1; ME(i).label = 'Idriss 2008';                  ME(i).func = @I2008;                      ME(i).mech=[0 0 1]; ME(i).type='regular';  ME(i).ref='https://doi.org/10.1193/1.2924362';
    i=i+1; ME(i).label = 'Chiou Youngs 2008';            ME(i).func = @CY2008;                     ME(i).mech=[0 0 1]; ME(i).type='regular';  ME(i).ref='https://doi.org/10.1193/1.2894832';
    i=i+1; ME(i).label = 'Boore Atkinson 2008';          ME(i).func = @BA2008;                     ME(i).mech=[0 0 1]; ME(i).type='regular';  ME(i).ref='https://doi.org/10.1193/1.2830434';
    i=i+1; ME(i).label = 'Campbell Bozorgnia 2008';      ME(i).func = @CB2008;                     ME(i).mech=[0 0 1]; ME(i).type='regular';  ME(i).ref='https://doi.org/10.1193/1.2857546';
    i=i+1; ME(i).label = 'Abrahamson Silva 2008';        ME(i).func = @AS2008;                     ME(i).mech=[0 0 1]; ME(i).type='regular';  ME(i).ref='https://doi.org/10.1193/1.2924360';
    i=i+1; ME(i).label = 'Abrahamson Silva 1997 (Horz)'; ME(i).func = @AS1997h;                    ME(i).mech=[0 0 1]; ME(i).type='regular';  ME(i).ref='https://doi.org/10.1785/gssrl.68.1.94';
    i=i+1; ME(i).label = 'Idriss 2014';                  ME(i).func = @I2014;                      ME(i).mech=[0 0 1]; ME(i).type='regular';  ME(i).ref ='https://doi.org/10.1193/070613EQS195M';
    i=i+1; ME(i).label = 'CY 2014';                      ME(i).func = @CY2014;                     ME(i).mech=[0 0 1]; ME(i).type='regular';  ME(i).ref='https://doi.org/10.1193/072813eqs219m';
    i=i+1; ME(i).label = 'CB 2014';                      ME(i).func = @CB2014;                     ME(i).mech=[0 0 1]; ME(i).type='regular';  ME(i).ref='https://doi.org/10.1193/062913eqs175m';
    i=i+1; ME(i).label = 'BSSA 2014';                    ME(i).func = @BSSA2014;                   ME(i).mech=[0 0 1]; ME(i).type='regular';  ME(i).ref='https://doi.org/10.1193/070113eqs184m';
    i=i+1; ME(i).label = 'ASK 2014';                     ME(i).func = @ASK2014;                    ME(i).mech=[0 0 1]; ME(i).type='regular';  ME(i).ref='https://doi.org/10.1193/070913eqs198m';
    i=i+1; ME(i).label = 'Akkar & Boomer 2007';          ME(i).func = @AkkarBoomer2007;            ME(i).mech=[0 0 1]; ME(i).type='regular';  ME(i).ref='https://doi.org/10.1002/eqe.679';
    i=i+1; ME(i).label = 'Akkar & Boomer 2010';          ME(i).func = @AkkarBoomer2010;            ME(i).mech=[0 0 1]; ME(i).type='regular';  ME(i).ref='https://doi.org/10.1785/gssrl.81.2.195';
    i=i+1; ME(i).label = 'Arroyo et al. 2010';           ME(i).func = @Arroyo2010;                 ME(i).mech=[1 0 0]; ME(i).type='regular';  ME(i).ref='https://doi.org/10.1007/s10950-010-9200-0';
    i=i+1; ME(i).label = 'Bindi et al. 2011';            ME(i).func = @Bindi2011;                  ME(i).mech=[0 0 1]; ME(i).type='regular';  ME(i).ref='https://doi.org/10.1007/s10518-011-9313-z';
    i=i+1; ME(i).label = 'Kanno et al. 2006';            ME(i).func = @Kanno2006;                  ME(i).mech=[1 1 0]; ME(i).type='regular';  ME(i).ref='http://dx.doi.org/10.1785/0120050138';
    i=i+1; ME(i).label = 'Cauzzi et al. 2015';           ME(i).func = @Cauzzi2015;                 ME(i).mech=[0 0 1]; ME(i).type='regular';  ME(i).ref='https://doi.org/10.1007/s10518-014-9685-y';
    i=i+1; ME(i).label = 'Du & Wang, 2012';              ME(i).func = @DW12;                       ME(i).mech=[0 0 1]; ME(i).type='regular';  ME(i).ref='https://doi.org/10.1002/eqe.2266';
    i=i+1; ME(i).label = 'Foulser-Piggott, Goda 2015';   ME(i).func = @FG15;                       ME(i).mech=[0 0 1]; ME(i).type='regular';  ME(i).ref='https://doi.org/10.1785/0120140316';
    i=i+1; ME(i).label = 'Travasarou et al. 2003';       ME(i).func = @TBA03;                      ME(i).mech=[0 0 1]; ME(i).type='regular';  ME(i).ref='https://doi.org/10.1002/eqe.270';
    i=i+1; ME(i).label = 'Bullock et al, 2017';          ME(i).func = @BU17;                       ME(i).mech=[1 1 1]; ME(i).type='regular';  ME(i).ref='https://doi.org/10.1785/0120160388';
    i=i+1; ME(i).label = 'Campbell,Bozorgnia 2010';      ME(i).func = @CB10;                       ME(i).mech=[0 0 1]; ME(i).type='regular';  ME(i).ref='https://doi.org/10.1193/1.3457158';
    i=i+1; ME(i).label = 'Campbell,Bozorgnia 2011';      ME(i).func = @CB11;                       ME(i).mech=[0 0 1]; ME(i).type='regular';  ME(i).ref='https://doi.org/10.1016/j.nucengdes.2011.04.020';
    i=i+1; ME(i).label = 'Campbell,Bozorgnia 2019';      ME(i).func = @CB19;                       ME(i).mech=[0 0 1]; ME(i).type='regular';  ME(i).ref='https://doi.org/10.1193/090818EQS212M';
    i=i+1; ME(i).label = 'Kramer & Mitchell, 2006';      ME(i).func = @KM06;                       ME(i).mech=[0 0 1]; ME(i).type='regular';  ME(i).ref='https://doi.org/10.1193/1.2194970';
    i=i+1; ME(i).label = 'PCE BCHydro (median)';         ME(i).func = @medianPCE_bchydro;          ME(i).mech=[1 1 0]; ME(i).type='regular';  ME(i).ref='www.google.com';
    
    i=i+1; ME(i).label = 'User Defined Model';           ME(i).func = @udm;                        ME(i).mech=[0 0 0]; ME(i).type='udm';      ME(i).ref='www.google.com';
    i=i+1; ME(i).label = 'Macedo et al. 2019';           ME(i).func = @Macedo2019;                 ME(i).mech=[1 1 0]; ME(i).type='cond';     ME(i).ref='https://doi.org/10.1785/0120180297';
    i=i+1; ME(i).label = 'Macedo et al. 2020';           ME(i).func = @Macedo2020;                 ME(i).mech=[0 0 1]; ME(i).type='cond';     ME(i).ref='www.google.com';
    i=i+1; ME(i).label = 'PCE NGA';                      ME(i).func = @PCE_nga;                    ME(i).mech=[0 0 1]; ME(i).type='pce';      ME(i).ref='www.google.com';
    i=i+1; ME(i).label = 'PCE BCHydro';                  ME(i).func = @PCE_bchydro;                ME(i).mech=[1 1 0]; ME(i).type='pce';      ME(i).ref='www.google.com';
    i=i+1; ME(i).label = 'Franky';                       ME(i).func = @franky;                     ME(i).mech=[0 0 0]; ME(i).type='frn';      ME(i).ref='www.google.com';

end

% magnitude scaling models
if group == 2
    i = 0;
    ME(1:4,1)=struct('label',[],'str',[],'func',[],'ref',[]);
    i=i+1;ME(i).label = 'Delta';                   ME(i).func = @delta;
    i=i+1;ME(i).label = 'Truncated Exponential';   ME(i).func = @truncexp;
    i=i+1;ME(i).label = 'Truncated Normal';        ME(i).func = @truncnorm;
    i=i+1;ME(i).label = 'Characteristic';          ME(i).func = @yc1985;
    i=i+1;ME(i).label = 'Table';                   ME(i).func = @magtable;
end

% spatial correlation models
if group ==3
    ME(1:4,1)=struct('label',[],'str',[],'func',[],'ref',[]);
    i = 0;
    i=i+1;ME(i).label = 'none';                              ME(i).func= @none_spatial;     ME(i).ref='www.google.com';
    i=i+1;ME(i).label = 'Jayaram N. and Baker J.W. (2009)';  ME(i).func= @JB_spatial_2009;  ME(i).ref='https://doi.org/10.1002/eqe.922';
    i=i+1;ME(i).label = 'Loth, C., and Baker, J. W. (2013)'; ME(i).func= @LB_spatial_2013;  ME(i).ref='https://doi.org/10.1002/eqe.2212';
    i=i+1;ME(i).label = 'Candia et al. (2019)';              ME(i).func= @SR_spatial_2019;  ME(i).ref='www.google.com';
end

% interperiod correlation models
if group == 4
    ME(1:13,1)=struct('label',[],'str',[],'func',[],'dependency',[],'ref',[]);                      % mechanism - magnitude - direction
    i = 0;
    i=i+1; ME(i).label = 'none';                     ME(i).func = @corr_none;                 ME(i).dependency = [0 0 0];    ME(i).ref='www.google.com';
    i=i+1; ME(i).label = 'Baker & Cornell 2006';     ME(i).func = @corr_BakerCornell2006;     ME(i).dependency = [0 0 1];    ME(i).ref='https://doi.org/10.1785/0120050060';
    i=i+1; ME(i).label = 'Baker & Jayaram 2008';     ME(i).func = @corr_BakerJayaram2008;     ME(i).dependency = [0 0 0];    ME(i).ref='https://doi.org/10.1193/1.2857544';
    i=i+1; ME(i).label = 'Jayaram et al. 2011';      ME(i).func = @corr_JayaramBaker2011;     ME(i).dependency = [1 0 0];    ME(i).ref='https://doi.org/10.12989/eas.2011.2.4.357';
    i=i+1; ME(i).label = 'Cimellaro 2013';           ME(i).func = @corr_Cimellaro2013;        ME(i).dependency = [0 0 1];    ME(i).ref='https://doi.org/10.1002/eqe.2248';
    i=i+1; ME(i).label = 'ASK2014 - NGA West2';      ME(i).func = @corr_Abrahamson2014;       ME(i).dependency = [0 0 0];    ME(i).ref='https://doi.org/10.1193/070913eqs198m';
    i=i+1; ME(i).label = 'Abrahamanson et al. 2016'; ME(i).func = @corr_BCHhydro2016;         ME(i).dependency = [0 0 0];    ME(i).ref='https://doi.org/10.1193/051712EQS188MR';
    i=i+1; ME(i).label = 'Baker & Bradley 2017';     ME(i).func = @corr_BakerBradley2017;     ME(i).dependency = [0 0 0];    ME(i).ref='https://doi.org/10.1193/060716EQS095M';
    i=i+1; ME(i).label = 'Jaimes & Candia 2019';     ME(i).func = @corr_JaimesCandia2019;     ME(i).dependency = [0 0 0];    ME(i).ref='www.google.com';
    i=i+1; ME(i).label = 'Candia et al. 2019';       ME(i).func = @corr_Candia2019;           ME(i).dependency = [1 0 0];    ME(i).ref='www.google.com';
    i=i+1; ME(i).label = 'Goda & Atkinson 2009';     ME(i).func = @corr_GodaAtkinson2009;     ME(i).dependency = [0 0 0];    ME(i).ref='https://doi.org/10.1785/0120090007';
    i=i+1; ME(i).label = 'Akkar & Sandikkaya 2014';  ME(i).func = @corr_Akkar2014;            ME(i).dependency = [0 0 0];    ME(i).ref='https://doi.org/10.1007/s10518-013-9537-1';
    i=i+1; ME(i).label = 'Ji et al. 2017';           ME(i).func = @corr_Ji2017;               ME(i).dependency = [0 1 0];    ME(i).ref='https://doi.org/10.1785/0120160291';
    i=i+1; ME(i).label = 'Campbell & Bozorgnia 2019';ME(i).func = @corr_CambpelBozorgnia2019; ME(i).dependency = [0 1 0];    ME(i).ref='https://doi.org/10.1193/090818EQS212M';
end

% psda
if group == 5
    ME(1:14,1)=struct('label',[],'str',[],'isregular',[],'func',[],'integrator',[],'primaryIM',[],'Safactor',[],'ref',[]);
    i=0;
    % subduction
    i=i+1;ME(i).label = 'BMT 2017 Sa(M)';      ME(i).func = @psda_BMT2017M;     ME(i).mechanism = 'subduction'; ME(i).integrator=1;  ME(i).primaryIM='SA(T=1.5Ts)';     ME(i).isregular=true;  ME(i).ref = 'https://doi.org/10.1061/(ASCE)GT.1943-5606.0001833';
    i=i+1;ME(i).label = 'BMT 2017 (PCE-M)';    ME(i).func = @psda_BMT2017_cdmM; ME(i).mechanism = 'subduction'; ME(i).integrator=6;  ME(i).primaryIM='SA(T=1.5Ts)';     ME(i).isregular=false; ME(i).ref = 'https://doi.org/10.1061/(ASCE)GT.1943-5606.0001833';
    
    % shallow crustal (bray et al)
    i=i+1;ME(i).label = 'BT 2007 Sa';          ME(i).func = @psda_BT2007;       ME(i).mechanism = 'crustal'; ME(i).integrator=2;  ME(i).primaryIM='SA(T=1.5Ts)';     ME(i).isregular=true;  ME(i).ref = 'https://doi.org/10.1061/(ASCE)1090-0241(2007)133:4(381)';
    i=i+1;ME(i).label = 'BT 2007 Sa(M)';       ME(i).func = @psda_BT2007M;      ME(i).mechanism = 'crustal'; ME(i).integrator=1;  ME(i).primaryIM='SA(T=1.5Ts)';     ME(i).isregular=true;  ME(i).ref = 'https://doi.org/10.1061/(ASCE)1090-0241(2007)133:4(381)';
    i=i+1;ME(i).label = 'BT 2007 (PCE)';       ME(i).func = @psda_BT2007_cdm;   ME(i).mechanism = 'crustal'; ME(i).integrator=5;  ME(i).primaryIM='SA(T=1.5Ts)';     ME(i).isregular=false; ME(i).ref = 'https://doi.org/10.1061/(ASCE)1090-0241(2007)133:4(381)';
    i=i+1;ME(i).label = 'BT 2007 (PCE-M)';     ME(i).func = @psda_BT2007_cdmM;  ME(i).mechanism = 'crustal'; ME(i).integrator=6;  ME(i).primaryIM='SA(T=1.5Ts)';     ME(i).isregular=false; ME(i).ref = 'https://doi.org/10.1061/(ASCE)1090-0241(2007)133:4(381)';
    i=i+1;ME(i).label = 'BM 2019 NonNF (M)';   ME(i).func = @psda_BM2019M;      ME(i).mechanism = 'crustal'; ME(i).integrator=1;  ME(i).primaryIM='SA(T=1.3Ts)';     ME(i).isregular=true;  ME(i).ref = 'https://www.ce.berkeley.edu/people/faculty/bray/research';
    %i=i+1;ME(i).label = 'BM 2019 NF';         ME(i).func = @psda_BM2019_NF;    ME(i).integrator=7;  ME(i).primaryIM='PGV-Sa(1.3Ts)'; ME(i).isregular=true;  ME(i).ref = 'http://www.google.com';
    
    % shallow crustal (other)
    i=i+1;ME(i).label = 'Jibson  2007 (M)';    ME(i).func = @psda_J07M;         ME(i).mechanism = 'crustal'; ME(i).integrator=1;  ME(i).primaryIM='PGA';           ME(i).isregular=true;  ME(i).ref = 'https://www.sciencedirect.com/science/article/pii/S0013795207000300?via%3Dihub';
    i=i+1;ME(i).label = 'Jibson  2007 Ia';     ME(i).func = @psda_J07Ia;        ME(i).mechanism = 'crustal'; ME(i).integrator=2;  ME(i).primaryIM='AI';            ME(i).isregular=true;  ME(i).ref = 'https://www.sciencedirect.com/science/article/pii/S0013795207000300?via%3Dihub';
    i=i+1;ME(i).label = 'RA 2011 (Rigid)';     ME(i).func = @psda_RA2011R;      ME(i).mechanism = 'crustal'; ME(i).integrator=3;  ME(i).primaryIM='PGV-PGA';       ME(i).isregular=true;  ME(i).ref = 'https://www.sciencedirect.com/science/article/pii/S0013795210002553';
    i=i+1;ME(i).label = 'RA 2011 (Flexible)';  ME(i).func = @psda_RA2011F;      ME(i).mechanism = 'crustal'; ME(i).integrator=4;  ME(i).primaryIM='PGV-PGA';       ME(i).isregular=true;  ME(i).ref = 'https://www.sciencedirect.com/science/article/pii/S0013795210002553';
    i=i+1;ME(i).label = 'RS 2009 (Scalar-M)';  ME(i).func = @psda_RS09M;        ME(i).mechanism = 'crustal'; ME(i).integrator=1;  ME(i).primaryIM='PGA';           ME(i).isregular=true;  ME(i).ref = 'http://www.nzsee.org.nz/db/Bulletin/Archive/42(1)0018.pdf';
    i=i+1;ME(i).label = 'RS 2009 (Vector)';    ME(i).func = @psda_RS09V;        ME(i).mechanism = 'crustal'; ME(i).integrator=3;  ME(i).primaryIM='PGV-PGA';       ME(i).isregular=true;  ME(i).ref = 'http://www.nzsee.org.nz/db/Bulletin/Archive/42(1)0018.pdf';
    i=i+1;ME(i).label = 'AM 1988';             ME(i).func = @psda_AM1988;       ME(i).mechanism = 'crustal'; ME(i).integrator=2;  ME(i).primaryIM='PGA';           ME(i).isregular=true;  ME(i).ref = 'https://doi.org/10.1002/eqe.4290160704';
    
    for j=1:length(ME)
        ME(j).Safactor=str2IM(regexp(ME(j).primaryIM,'\-','split'));
        ME(j).Safactor=ME(j).Safactor(:)';
    end
    
end

% settlement and tilt models
if group == 6
    ME(1:4,1)=struct('label',[],'func',[],'str',[],'IM',[],'integrator',[],'ref',[]);
    i=0;
    i=i+1;ME(i).label = 'Bullock 2018 (Settlement)';            ME(i).func = @set_B18;   ME(i).IM='CAV';             ME(i).integrator=1; ME(i).ref = 'https://doi.org/10.1680/jgeot.17.P.174';
    i=i+1;ME(i).label = 'Ishihara 2017 (Settlement)';           ME(i).func = @set_I17;   ME(i).IM='PGA-SA(T=1)-CAV'; ME(i).integrator=2; ME(i).ref = 'https://doi.org/10.1016/j.soildyn.2017.08.026';
    i=i+1;ME(i).label = 'Bullock, 2018 (Tilt Empirical)';       ME(i).func = @til_B18e;  ME(i).IM='CAV';             ME(i).integrator=1; ME(i).ref = 'https://doi.org/10.1680/jgeot.17.P.174';
    i=i+1;ME(i).label = 'Bullock, 2018 (Tilt Semi Empirical)';  ME(i).func = @til_B18s;  ME(i).IM='CAV-VGI';         ME(i).integrator=1; ME(i).ref = 'https://doi.org/10.1680/jgeot.17.P.174';
end


for i=1:length(ME)
    ME(i).str=func2str(ME(i).func);
end


%% selective method
if nargin==2
    val = varargin{1};
    ME = ME(val);
end