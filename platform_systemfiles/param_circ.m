function [param,rate]=param_circ(r0,source,ellipsoid,hparam)

RA      = source.numgeom(3:4);
gmm     = source.gmm;
M       = source.mscl(:,1);
nM      = size(M,1);
rf      = source.hypm;
nR      = size(rf,1);
media   = source.media;
[iR,iM] = meshgrid(1:nR,1:nM);
M       = M(iM(:));
rf      = rf(iR(:),:);

[~,indVs30]=intersect(hparam,'VS30'); VS30 = media(indVs30);

%% EVALUATES RUPTURE AREA AND SOURCE NORMAL VECTOR
rupArea   = rupRelation(M,0,RA(1));
if size(source.normal,1)==1
    n = zeros(length(iR(:)),3);
    n(:,1)=source.normal(1);
    n(:,2)=source.normal(2);
    n(:,3)=source.normal(3);
else
    n         = source.normal(iR(:),:);
end
rateM     = source.mscl(iM(:),2);
rateR     = source.aream(iR(:))/sum(source.aream);
rate      = rateM.*rateR/(rateM'*rateR);

% this is required if you use NGAWest2 models on non-rectangular sources
rrup  = dist_rrup(r0,rf,rupArea,n);
UNK   = ones(size(rrup))*999;
if gmm.Rmetric(2), rhyp  = dist_rhyp(r0,rf);           end
if gmm.Rmetric(3), rjb   = dist_rjb(r0,rf,rupArea,n,ellipsoid);  end
if gmm.Rmetric(8), zhyp  = dist_zhyp(r0,rf,ellipsoid); end
if gmm.Rmetric(9), ztor  = dist_ztor(r0,rf,rupArea,n,ellipsoid);end

%% GMM Parameters
Ndepend  = 1;
isfrn    = false;
switch gmm.type
    case 'regular', str_test = func2str(gmm.handle);
    case 'cond',    str_test = func2str(gmm.cond);
    case 'udm' ,    str_test = 'udm';
    case 'pce' ,    str_test = func2str(gmm.handle);
    case 'frn'
        isfrn   = true;
        Ndepend = length(gmm.usp);
        funcs   = cell(1,Ndepend);
        IMlist  = cell(1,Ndepend);
        PARAM   = cell(1,Ndepend);
end

for jj=1:Ndepend
    if isfrn
        if strcmp(gmm.usp{jj}.type,'regular')
            str_test = func2str(gmm.usp{jj}.handle);
        else
            str_test = func2str(gmm.usp{jj}.cond);
        end
        usp      = gmm.usp{jj}.usp;
    else
        usp     = gmm.usp;
    end
    
    dip = 90;
    switch str_test
        case 'Youngs1997',                    param = [M,rrup,zhyp,VS30,usp];
        case 'AtkinsonBoore2003',             param = [M,rrup,zhyp,VS30,usp];
        case 'Zhao2006',                      param = [M,rrup,zhyp,VS30,usp];
        case 'Mcverry2006',                   param = [M,rrup,zhyp,VS30,usp];
        case 'ContrerasBoroschek2012',        param = [M,rrup,zhyp,VS30,usp];
        case 'BCHydro2012',                   param = [M,rrup,rhyp,zhyp,VS30,usp];
        case 'BCHydro2018',                   param = [M,rrup,ztor,VS30,usp];
        case 'Kuehn2020',                     param = [M,rrup,ztor,VS30,usp];
        case 'Parker2020',                    param = [M,rrup,zhyp,VS30,usp];
        case 'Arteta2018',                    param = [M,rhyp,VS30,usp];
        case 'Idini2016',                     param = [M,rrup,rhyp,zhyp,VS30,usp];
        case 'MontalvaBastias2017',           param = [M,rrup,rhyp,zhyp,VS30,usp];
        case 'MontalvaBastias2017HQ',         param = [M,rrup,rhyp,zhyp,VS30,usp];
        case 'Montalva2018'
            [~,indf0] = intersect(hparam,'f0');    f0    = media(indf0);
            param = [M,rrup,rhyp,zhyp,VS30,f0,usp];
        case 'SiberRisk2019',                 param = [M,rrup,rhyp,zhyp,VS30,usp];
        case 'Garcia2005',                    param = [M,rrup,rhyp,zhyp,usp];
        case 'Jaimes2006',                    param = [M,rrup,usp];
        case 'Jaimes2015',                    param = [M,rrup,usp];
        case 'Jaimes2016',                    param = [M,rrup,usp];
        case 'GarciaJaimes2017',              param = [M,rrup,usp];
        case 'GarciaJaimes2017HV',            param = [M,rrup,usp];
        case 'Bernal2014',                    param = [M,rrup,zhyp,usp];
        case 'Sadigh1997',                    param = [M,rrup,VS30,usp];
        case 'I2008',                         param = [M,rrup,VS30,usp];
        case 'CY2008',                        param = [M,rrup,UNK,-UNK,ztor,dip,VS30,usp];
        case 'BA2008',                        param = [M,rjb,VS30,usp];
        case 'CB2008',                        param = [M,rrup,rjb,ztor, nan,VS30,usp];
        case 'AS2008',                        param = [M,rrup,UNK,-UNK,ztor,dip,999,VS30,usp];
        case 'AS1997h',                       param = [M,rrup,VS30,usp];
        case 'I2014',                         param = [M,rrup,VS30,usp];
        case 'CY2014',                        param = [M,rrup,UNK,UNK,ztor,dip,VS30,usp];
        case 'CB2014',                        param = [M,rrup,UNK,UNK,UNK,ztor,'unk',dip,999,VS30,usp];
        case 'BSSA2014',                      param = [M,rjb,VS30,usp];
        case 'ASK2014',                       param = [M,rrup,UNK,-UNK,UNK,ztor,dip,999,VS30,usp];
        case 'AkkarBoomer2007',               param = [M,rjb,usp];
        case 'AkkarBoomer2010',               param = [M,rjb,usp];
        case 'Arroyo2010',                    param = [M,rrup,usp];
        case 'Bindi2011',                     param = [M,rjb,VS30,usp];
        case 'Kanno2006',                     param = [M,rrup,zhyp,VS30,usp];
        case 'Cauzzi2015',                    param = [M,rrup,rhyp,VS30,usp];
        case 'DW12',                          param = [M,rrup,usp];
        case 'FG15',                          param = [M,rrup,zhyp,VS30,usp];
        case 'TBA03',                         param = [M,rrup,usp];
        case 'BU17',                          param = [M,rrup,zhyp,usp];
        case 'CB10',                          param = [M,rrup,rjb,ztor,dip,VS30,usp];
        case 'CB11',                          param = [M,rrup,rjb,ztor,dip,VS30,usp];
        case 'CB19',                          param = [M,rrup,rjb,UNK,UNK,ztor,zhyp,dip,999,VS30,usp];
        case 'KM06',                          param = [M,rrup,usp];
        case 'medianPCE_bchydro',             param = [M,rrup,VS30,usp];
        case 'PCE_nga',                       param = [M,rrup,VS30,usp];
        case 'PCE_bchydro',                   param = [M,rrup,VS30,usp];
        case 'udm'
            var      = gmm.var;
            txt      = regexp(var.syntax,'\(','split');
            args     = regexp(txt{2}(1:end-1),'\,','split');
            args     = strtrim(args);
            args(1)  = [];
            param    = cell(1,4+length(args));
            param{1} = str2func(strtrim(txt{1}));
            param{2} = var.vector;
            param{3} = var.residuals;
            uspcont  = 2;
            for cont=1:length(args)
                f = var.(args{cont});
                if strcmpi(f.tag,'magnitude')
                    param{4+cont}=M;
                end
                
                if strcmpi(f.tag,'distance')
                    fval = find(f.value);
                    switch fval
                        case 1 , param{4+cont}=rrup;
                        case 2 , param{4+cont}=rhyp;
                        case 3 , param{4+cont}=rjb;
                        case 4 , param{4+cont}=repi;
                        case 6 , param{4+cont}=rx;
                        case 7 , param{4+cont}=ry0;
                        case 8 , param{4+cont}=zhyp;
                        case 9 , param{4+cont}=rztor;
                        case 10, param{4+cont}=rzbor;
                        case 11, param{4+cont}=rzbot;
                    end
                end
                
                if strcmpi(f.tag,'Vs30')
                    param{4+cont} = source.media;
                    uspcont=uspcont+1;
                end
                
                if strcmpi(f.tag,'param')
                    switch f.type
                        case 'string'
                            param{4+cont}=usp{uspcont};
                        case 'double'
                            param{4+cont}=str2double(usp{uspcont});
                            
                    end
                    uspcont=uspcont+1;
                end
                
            end
    end
    
    switch gmm.type
        case 'cond'
            switch func2str(gmm.handle)
                case 'Macedo2019', param = {M,rrup,gmm.txt{4},gmm.txt{5},VS30,gmm.cond,param{:}}; %#ok<CCAT>
                case 'Macedo2020', param = {M,rrup,VS30         ,gmm.cond,param{:}}; %#ok<CCAT>
            end
        case 'frn'
            funcs {jj}=gmm.usp{jj}.handle;
            IMlist{jj}=gmm.usp{jj}.T;
            
            switch gmm.usp{jj}.type
                case 'regular'
                    PARAM {jj}=param;
                case 'cond'
                    switch func2str(gmm.usp{jj}.handle)
                        case 'Macedo2019', PARAM{jj} = {M,rrup,gmm.usp{jj}.txt{4},gmm.usp{jj}.txt{5},VS30,gmm.usp{jj}.cond,param{:}}; %#ok<CCAT>
                        case 'Macedo2020', PARAM{jj} = {M,rrup,VS30,gmm.usp{jj}.cond,param{:}}; %#ok<CCAT>
                    end
            end
    end
end

if isfrn
    param=[funcs,IMlist,PARAM];
end

