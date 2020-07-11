function[val,units]=isIMadmisible(To,gmm,rangeSA,rangeSV,rangeSD,rangeHV)

switch gmm
    case 'Youngs1997'                , IMlist=[1 0 0 0 0 0 0 1 0 0 0];    cf=[1,1,1,1,1,1,1,1,1,1,1];
    case 'AtkinsonBoore2003'         , IMlist=[1 0 0 0 0 0 0 1 0 0 0];    cf=[1/980.66,1,1,1,1,1,1,1/980.66,1,1,1];
    case 'Zhao2006'                  , IMlist=[1 0 0 0 0 0 0 1 0 0 0];    cf=[1/980.66,1,1,1,1,1,1,1/980.66,1,1,1];
    case 'Mcverry2006'               , IMlist=[1 0 0 0 0 0 0 1 0 0 0];    cf=[1,1,1,1,1,1,1,1,1,1,1];
    case 'ContrerasBoroschek2012'    , IMlist=[1 0 0 0 0 0 0 1 0 0 0];    cf=[1,1,1,1,1,1,1,1,1,1,1];
    case 'BCHydro2012'               , IMlist=[1 0 0 0 0 0 0 1 0 0 0];    cf=[1,1,1,1,1,1,1,1,1,1,1];
    case 'BCHydro2018'               , IMlist=[1 0 0 0 0 0 0 1 0 0 0];    cf=[1,1,1,1,1,1,1,1,1,1,1];
    case 'Kuehn2020'                 , IMlist=[1 1 0 0 0 0 0 1 0 0 0];    cf=[1,1,1,1,1,1,1,1,1,1,1];
    case 'Parker2020'                , IMlist=[1 1 0 0 0 0 0 1 0 0 0];    cf=[1,1,1,1,1,1,1,1,1,1,1];
    case 'Arteta2018'                , IMlist=[1 0 0 0 0 0 0 1 0 0 0];    cf=[1,1,1,1,1,1,1,1,1,1,1];
    case 'Idini2016'                 , IMlist=[1 0 0 0 0 0 0 1 0 0 0];    cf=[1,1,1,1,1,1,1,1,1,1,1];
    case 'MontalvaBastias2017'       , IMlist=[1 0 0 0 0 0 0 1 0 0 0];    cf=[1,1,1,1,1,1,1,1,1,1,1];
    case 'MontalvaBastias2017HQ'     , IMlist=[1 0 0 0 0 0 0 1 0 0 0];    cf=[1,1,1,1,1,1,1,1,1,1,1];
    case 'Montalva2018'              , IMlist=[1 0 0 0 0 0 0 1 0 0 0];    cf=[1,1,1,1,1,1,1,1,1,1,1];    
    case 'SiberRisk2019'             , IMlist=[1 1 0 0 0 0 0 1 0 0 0];    cf=[1,100,1,1,1,1,1,1,1,1,1];
    case 'Garcia2005'                , IMlist=[1 1 0 0 0 0 0 1 0 0 0];    cf=[1/980.66,1,1,1,1,1,1,1/980.66,1,1,1];
    case 'Jaimes2006'                , IMlist=[1 0 0 0 0 0 0 1 0 0 0];    cf=[1/980.66,1,1,1,1,1,1,1/980.66,1,1,1];
    case 'Jaimes2015'                , IMlist=[0 1 0 0 0 0 0 1 0 0 0];    cf=[1/980.66,1,1,1,1,1,1,1/980.66,1,1,1];
    case 'Jaimes2016'                , IMlist=[1 1 0 0 0 0 0 1 0 0 0];    cf=[1/980.66,1,1,1,1,1,1,1/980.66,1,1,1];
    case 'GarciaJaimes2017'          , IMlist=[1 1 0 0 0 0 0 1 0 0 0];    cf=[1/980.66,1,1,1,1,1,1,1/980.66,1,1,1];
    case 'GarciaJaimes2017HV'        , IMlist=[0 0 0 0 0 0 0 0 0 0 1];    cf=[1,1,1,1,1,1,1,1,1,1,1];
    case 'Bernal2014'                , IMlist=[1 0 0 0 0 0 0 1 0 0 0];    cf=[1/980.66,1,1,1,1,1,1,1/980.66,1,1,1];
    case 'Sadigh1997'                , IMlist=[1 0 0 0 0 0 0 1 0 0 0];    cf=[1,1,1,1,1,1,1,1,1,1,1];
    case 'I2008'                     , IMlist=[1 0 0 0 0 0 0 1 0 0 0];    cf=[1,1,1,1,1,1,1,1,1,1,1];
    case 'CY2008'                    , IMlist=[1 1 0 0 0 0 0 1 0 0 0];    cf=[1,1,1,1,1,1,1,1,1,1,1];
    case 'BA2008'                    , IMlist=[1 1 0 0 0 0 0 1 0 0 0];    cf=[1,1,1,1,1,1,1,1,1,1,1];
    case 'CB2008'                    , IMlist=[1 1 1 0 0 0 0 1 0 0 0];    cf=[1,1,1,1,1,1,1,1,1,1,1];
    case 'AS2008'                    , IMlist=[1 1 0 0 0 0 0 1 0 0 0];    cf=[1,1,1,1,1,1,1,1,1,1,1];
    case 'AS1997h'                   , IMlist=[0 0 0 0 0 0 0 1 0 0 0];    cf=[1,1,1,1,1,1,1,1,1,1,1];
    case 'I2014'                     , IMlist=[1 0 0 0 0 0 0 1 0 0 0];    cf=[1,1,1,1,1,1,1,1,1,1,1];
    case 'CY2014'                    , IMlist=[1 1 0 0 0 0 0 1 0 0 0];    cf=[1,1,1,1,1,1,1,1,1,1,1];
    case 'CB2014'                    , IMlist=[1 1 0 0 0 0 0 1 0 0 0];    cf=[1,1,1,1,1,1,1,1,1,1,1];
    case 'BSSA2014'                  , IMlist=[1 1 0 0 0 0 0 1 0 0 0];    cf=[1,1,1,1,1,1,1,1,1,1,1];
    case 'ASK2014'                   , IMlist=[1 1 0 0 0 0 0 1 0 0 0];    cf=[1,1,1,1,1,1,1,1,1,1,1];
    case 'AkkarBoomer2007'           , IMlist=[1 0 0 0 0 0 0 0 0 1 0];    cf=[1/980.66,1,1,1,1,1,1,1,1,1,1];
    case 'AkkarBoomer2010'           , IMlist=[1 1 0 0 0 0 0 1 0 0 0];    cf=[1/980.66,1,1,1,1,1,1,1/980.66,1,1,1];
    case 'Arroyo2010'                , IMlist=[1 0 0 0 0 0 0 1 0 0 0];    cf=[1/980.66,1,1,1,1,1,1,1/980.66,1,1,1];
    case 'Bindi2011'                 , IMlist=[1 1 0 0 0 0 0 1 0 0 0];    cf=[1/980.66,1,1,1,1,1,1,1/980.66,1,1,1];
    case 'Kanno2006'                 , IMlist=[1 1 0 0 0 0 0 1 0 0 0];    cf=[1/980.66,1,1,1,1,1,1,1/980.66,1,1,1];
    case 'Cauzzi2015'                , IMlist=[1 1 0 0 0 0 0 1 0 1 0];    cf=[1/980.66,1,1,1,1,1,1,1/980.66,1,1,1];
    case 'DW12'                      , IMlist=[0 0 0 0 1 0 0 0 0 0 0];    cf=[1,1,1,1,980.66,1,1,1,1,1,1];
    case 'FG15'                      , IMlist=[0 0 0 0 1 1 0 0 0 0 0];    cf=[1,1,1,1,100,100,1,1,1,1,1];
    case 'TBA03'                     , IMlist=[0 0 0 0 0 1 0 0 0 0 0];    cf=[1,1,1,1,1,100,1,1,1,1,1];
    case 'BU17'                      , IMlist=[0 0 0 0 1 1 1 0 0 0 0];    cf=[1,1,1,1,1,1,1,1,1,1,1];
    case 'CB10'                      , IMlist=[0 0 0 0 1 0 0 0 0 0 0];    cf=[1,1,1,1,980.66,1,1,1,1,1,1];
    case 'CB11'                      , IMlist=[0 0 0 0 1 0 0 0 0 0 0];    cf=[1,1,1,1,980.66,1,1,1,1,1,1];
    case 'CB19'                      , IMlist=[0 0 0 0 1 1 0 0 0 0 0];    cf=[1,1,1,1,100,100,1,1,1,1,1];
    case 'KM06'                      , IMlist=[0 0 0 0 1 0 0 0 0 0 0];    cf=[1,1,1,1,100,1,1,1,1,1,1];
    case 'medianPCE_bchydro'         , IMlist=[1 0 0 0 1 0 0 0 0 0 0];    cf=[1,1,1,1,100,1,1,1,1,1,1];
    case 'Macedo2019'                , IMlist=[0 0 0 0 0 1 0 0 0 0 0];    cf=[1 1 1 1 1 100 1 1 1 1 1];
    case 'Macedo2020'                , IMlist=[1 0 0 0 0 0 0 1 0 0 0];    cf=[1,1,1,1,1,1,1,1,1,1,1];
end

val=false;
rTo=real(To);
if imag(To)==0 && To<=0
    % Non-spectral GMMs
    switch To
        case  0, val = IMlist(1)==1; units = cf(1);
        case -1, val = IMlist(2)==1; units = cf(2);
        case -2, val = IMlist(3)==1; units = cf(3);
        case -3, val = IMlist(4)==1; units = cf(4);
        case -4, val = IMlist(5)==1; units = cf(5);
        case -5, val = IMlist(6)==1; units = cf(6);
        case -6, val = IMlist(7)==1; units = cf(7);
    end
    % spectral GMMs
elseif imag(To)==0 && To>0,  val = all([IMlist( 8),rTo>=rangeSA(1),rTo<=rangeSA(end)]); units = cf(8);
elseif imag(To)==1 && rTo>0, val = all([IMlist( 9),rTo>=rangeSV(1),rTo<=rangeSV(end)]); units = cf(9);
elseif imag(To)==2 && rTo>0, val = all([IMlist(10),rTo>=rangeSD(1),rTo<=rangeSD(end)]); units = cf(10);
elseif imag(To)==3 && rTo>0, val = all([IMlist(11),rTo>=rangeHV(1),rTo<=rangeHV(end)]); units = cf(11);
end

