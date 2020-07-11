function[Rrup,param]=mGMPErrupLoop(fun,param,SUB,SC)

str  = func2str(fun);
switch str
    case 'Youngs1997'
        param{1} = SUB.rpMag;
        param{2} = SUB.rpRrup;
        param{3} = SUB.rpZhyp;
        Rrup     = SUB.rpRrup;
        
    case 'AtkinsonBoore2003'
        param{1} = SUB.rpMag;
        param{2} = SUB.rpRrup;
        param{3} = SUB.rpZhyp;
        Rrup     = SUB.rpRrup;
        
    case 'Zhao2006'
        param{1} = SUB.rpMag;
        param{2} = SUB.rpRrup;
        param{3} = SUB.rpZhyp;
        Rrup     = SUB.rpRrup;
        
    case 'Mcverry2006'
        param{1} = SUB.rpMag;
        param{2} = SUB.rpRrup;
        param{3} = SUB.rpZhyp;
        Rrup     = SUB.rpRrup;
        
    case 'ContrerasBoroschek2012'
        param{1} = SUB.rpMag;
        param{2} = SUB.rpRrup;
        param{3} = SUB.rpZhyp;
        Rrup     = SUB.rpRrup;
        
    case 'BCHydro2012'
        param{1} = SUB.rpMag;
        param{2} = SUB.rpRrup;
        param{3} = SUB.rpRhyp;
        param{4} = SUB.rpZhyp;
        Rrup     = SUB.rpRrup;
        
    case 'BCHydro2018'
        param{1} = SUB.rpMag;
        param{2} = SUB.rpRrup;
        param{3} = SUB.rpZtor;
        Rrup     = SUB.rpRrup;
        
    case 'Kuehn2020'
        param{1} = SUB.rpMag;
        param{2} = SUB.rpRrup;
        param{3} = SUB.rpZtor;
        Rrup     = SUB.rpRrup;
        
    case 'Parker2020'
        param{1} = SUB.rpMag;
        param{2} = SUB.rpRrup;
        param{3} = SUB.rpZhyp;
        Rrup     = SUB.rpRrup;
        
    case 'Arteta2018'
        param{1} = SUB.rpMag;
        param{2} = SUB.rpRrup;
        Rrup     = SUB.rpRrup;
        
    case 'Idini2016'
        param{1} = SUB.rpMag;
        param{2} = SUB.rpRrup;
        param{3} = SUB.rpRhyp;
        param{4} = SUB.rpZhyp;
        Rrup     = SUB.rpRrup;
        
    case 'MontalvaBastias2017'
        param{1} = SUB.rpMag;
        param{2} = SUB.rpRrup;
        param{3} = SUB.rpRhyp;
        param{4} = SUB.rpZhyp;
        Rrup     = SUB.rpRrup;
        
    case 'MontalvaBastias2017HQ'
        param{1} = SUB.rpMag;
        param{2} = SUB.rpRrup;
        param{3} = SUB.rpRhyp;
        param{4} = SUB.rpZhyp;
        Rrup     = SUB.rpRrup;
        
    case 'Montalva2018'
        param{1} = SUB.rpMag;
        param{2} = SUB.rpRrup;
        param{3} = SUB.rpRhyp;
        param{4} = SUB.rpZhyp;
        Rrup     = SUB.rpRrup;
        
    case 'SiberRisk2019'
        param{1} = SUB.rpMag;
        param{2} = SUB.rpRrup;
        param{3} = SUB.rpRhyp;
        param{4} = SUB.rpZhyp;
        Rrup     = SUB.rpRrup;
        
    case 'Garcia2005'
        param{1} = SUB.rpMag;
        param{2} = SUB.rpRrup;
        param{3} = SUB.rpRhyp;
        param{4} = SUB.rpZhyp;
        Rrup     = SUB.rpRrup;
        
    case 'Jaimes2006'
        param{1} = SUB.rpMag;
        param{2} = SUB.rpRrup;
        Rrup     = SUB.rpRrup;
        
    case 'Jaimes2015'
        param{1} = SUB.rpMag;
        param{2} = SUB.rpRrup;
        Rrup     = SUB.rpRrup;
        
    case 'Jaimes2016'
        param{1} = SUB.rpMag;
        param{2} = SUB.rpRrup;
        Rrup     = SUB.rpRrup;
        
    case 'GarciaJaimes2017'
        param{1} = SUB.rpMag;
        param{2} = SUB.rpRrup;
        Rrup     = SUB.rpRrup;
        
    case 'GarciaJaimes2017HV'
        param{1} = SUB.rpMag;
        param{2} = SUB.rpRrup;
        Rrup     = SUB.rpRrup;
        
    case 'Bernal2014'
        param{1} = SUB.rpMag;
        param{2} = SUB.rpRrup;
        param{3} = SUB.rpZhyp;
        Rrup     = SUB.rpRrup;
        
    case 'Sadigh1997'
        param{1} = SC.rpMag;
        param{2} = SC.rpRrup;
        Rrup     = SC.rpRrup;
        
    case 'I2008'
        param{1} = SC.rpMag;
        param{2} = SC.rpRrup;
        Rrup     = SC.rpRrup;
        
    case 'CY2008'
        param{1} = SC.rpMag;
        param{2} = SC.rpRrup;
        param{3} = SC.rpRjb;
        param{4} = SC.rpRx;
        param{5} = SC.rpZtor;
        Rrup     = SC.rpRrup;
        
    case 'BA2008'
        param{1} = SC.rpMag;
        param{2} = SC.rpRjb;
        Rrup     = SC.rpRrup;
        
    case 'CB2008'
        param{1} = SC.rpMag;
        param{2} = SC.rpRrup;
        param{3} = SC.rpRjb;
        param{4} = SC.rpZtor;
        Rrup     = SC.rpRrup;
        
    case 'AS2008'
        param{1} = SC.rpMag;
        param{2} = SC.rpRrup;
        param{3} = SC.rpRjb;
        param{4} = SC.rpRx;
        param{5} = SC.rpZtor;
        Rrup     = SC.rpRrup;
        
    case 'AS1997h'
        param{1} = SC.rpMag;
        param{2} = SC.rpRrup;
        Rrup     = SC.rpRrup;
        
    case 'I2014'
        param{1} = SC.rpMag;
        param{2} = SC.rpRrup;
        Rrup     = SC.rpRrup;
        
    case 'CY2014'
        param{1} = SC.rpMag;
        param{2} = SC.rpRrup;
        param{3} = SC.rpRjb;
        param{4} = SC.rpRx;
        param{5} = SC.rpZtor;
        Rrup     = SC.rpRrup;
        
    case 'CB2014'
        param{1} = SC.rpMag;
        param{2} = SC.rpRrup;
        param{3} = SC.rpRjb;
        param{4} = SC.rpRx;
        param{5} = SC.rpZhyp;
        param{6} = SC.rpZtor;
        param{7} = SC.Zbot;
        Rrup     = SC.rpRrup;
        
    case 'BSSA2014'
        param{1} = SC.rpMag;
        param{2} = SC.rpRjb;
        Rrup     = SC.rpRrup;
        
    case 'ASK2014'
        param{1} = SC.rpMag;
        param{2} = SC.rpRrup;
        param{3} = SC.rpRjb;
        param{4} = SC.rpRx;
        param{5} = SC.rpRy0;
        param{6} = SC.rpZtor;
        Rrup     = SC.rpRrup;
        
    case 'AkkarBoomer2007'
        param{1} = SC.rpMag;
        param{2} = SC.rpRjb;
        Rrup     = SC.rpRrup;
        
    case 'AkkarBoomer2010'
        param{1} = SC.rpMag;
        param{2} = SC.rpRjb;
        Rrup     = SC.rpRrup;
        
    case 'Arroyo2010'
        param{1} = SC.rpMag;
        param{2} = SC.rpRrup;
        Rrup     = SC.rpRrup;
        
    case 'Bindi2011'
        param{1} = SC.rpMag;
        param{2} = SC.rpRjb;
        Rrup     = SC.rpRrup;
        
    case 'Kanno2006'
        param{1} = SC.rpMag;
        param{2} = SC.rpRrup;
        param{3} = SC.rpZhyp;
        Rrup     = SC.rpRrup;
        
    case 'Cauzzi2015'
        param{1} = SC.rpMag;
        param{2} = SC.rpRrup;
        param{3} = SC.rpRhyp;
        Rrup     = SC.rpRrup;
        
    case 'DW12'
        param{1} = SC.rpMag;
        param{2} = SC.rpRrup;
        Rrup     = SC.rpRrup;
        
    case 'FG15'
        param{1} = SC.rpMag;
        param{2} = SC.rpRrup;
        param{3} = SC.rpZhyp;
        Rrup     = SC.rpRrup;
        
    case 'TBA03'
        param{1} = SC.rpMag;
        param{2} = SC.rpRrup;
        Rrup     = SC.rpRrup;
        
    case 'BU17'
        param{1} = SC.rpMag;
        param{2} = SC.rpRrup;
        param{3} = SC.rpZhyp;
        Rrup     = SC.rpRrup;
        
    case 'CB10'
        param{1} = SC.rpMag;
        param{2} = SC.rpRrup;
        param{3} = SC.rpRjb;
        param{4} = SC.rpZtor;
        Rrup     = SC.rpRrup;
        
    case 'CB11'
        param{1} = SC.rpMag;
        param{2} = SC.rpRrup;
        param{3} = SC.rpRjb;
        param{4} = SC.rpZtor;
        Rrup     = SC.rpRrup;
        
    case 'CB19'
        param{1} = SC.rpMag;
        param{2} = SC.rpRrup;
        param{3} = SC.rpRjb;
        param{4} = SC.rpRx;
        param{5} = SC.rpZhyp;
        param{6} = SC.rpZtor;
        param{7} = SC.Zbot;
        Rrup     = SC.rpRrup;
        
    case 'KM06'
        param{1} = SC.rpMag;
        param{2} = SC.rpRrup;
        Rrup     = SC.rpRrup;
        
    case 'medianPCE_bchydro'
        param{1} = SUB.rpMag;
        param{2} = SUB.rpRrup;
        Rrup     = SUB.rpRrup;        
        
end