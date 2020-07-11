function [param,rate]=param_rect(r0,source,ellipsoid,hparam)

if isempty(ellipsoid.Code)
    r0=r0([2 1 3]);
end

Area    = source.aream;
MM      = source.mscl(:,1);
rateM   = source.mscl(:,2);
gmm     = source.gmm;
RA      = source.numgeom(3:4);
L       = sum(abs(source.p(:,1)))/2;%source.numgeom(8);
dip     = source.numgeom(7);
if source.numgeom(9)==-999
    width = Area/L;
else
    width   = source.numgeom(9);
end

aratio  = source.numgeom(10);
spacing = source.numgeom(11);
Rmetric = gmm.Rmetric;
media   = source.media;
[~,indVs30]=intersect(hparam,'VS30'); VS30 = media(indVs30);

%% RUPTURE AREA AND SCENARIOS
if RA(2)==0
    rupArea   = rupRelation(MM,0,RA(1));
else
    NA = 25;
    x  = linspace(-2,2,NA)';
    dx = x(2)-x(1);
    rateA = normcdf(x+dx/2,0,1)-normcdf(x-dx/2,0,1);
    rateA = rateA/sum(rateA);
    rupArea = zeros(length(MM),NA);
    for i=1:NA
        rupArea(:,i) = rupRelation(MM,x(i),RA(1));
    end
    rupArea = min(rupArea,Area);
    [iM,iA]=meshgrid(1:length(MM),1:length(x));
    MM     = MM(iM(:));
    rateM = rateM(iM(:)).*rateA(iA(:));
    rupArea = rupArea';
    rupArea = rupArea(iA(:));
end

rupWidth  = min(sqrt(rupArea/aratio),width);  % preserve area at expense of aspect ratio
rupLength = min(rupArea./rupWidth,L);     % preserve area at expense of aspect ratio

nM    = length(MM);
if Rmetric(1),  rrup  = cell(1,nM);end
if Rmetric(2),  rhyp  = cell(1,nM);end
if Rmetric(3),  rjb   = cell(1,nM);end
if Rmetric(6),  rx    = cell(1,nM);end
if Rmetric(7),  ry0   = cell(1,nM);end
if Rmetric(8),  zhyp  = cell(1,nM);end
if Rmetric(9),  ztor  = cell(1,nM);end

p     = source.p;
pmean = source.hypm;
rot   = source.normal;

geom.pmean  = pmean;
geom.normal =-rot(:,3)';
geom.rot    = rot;
geom.dip    = dip;

for i=1:nM
    RL    = rupLength(i);
    RW    = rupWidth (i);
    xmin  = p(1,1)+RL/2;
    xmax  = p(2,1)-RL/2;
    ymin  = p(2,2)+RW/2;
    ymax  = p(3,2)-RW/2;
    if xmin>xmax
        xavg = 1/2*(xmin+xmax);
        xmin = xavg;
        xmax = xavg;
    end
    
    if ymin>ymax
        yavg = 1/2*(ymin+ymax);
        ymin = yavg;
        ymax = yavg;
    end
    
    dx    = max(xmax-xmin,0);
    dy    = max(ymax-ymin,0);
    NX    = ceil(dx/spacing)+1;
    NY    = ceil(dy/spacing)+1;
    locx  = linspace(xmin,xmax,NX);
    locy  = linspace(ymin,ymax,NY);
    [locx,locy] = meshgrid(locx,locy);
    nR    = numel(locx);
    locxy = [locx(:),locy(:),zeros(nR,1)]; % Y,X,Z coordinates of
    rf      = bsxfun(@plus,pmean,locxy*rot');
    
    if Rmetric(1),  rrup{i}  = dist_rrup4  (r0,rf,RW,RL,geom);end
    if Rmetric(2),  rhyp{i}  = dist_rhyp4  (r0,rf,RW,RL,geom);end
    if Rmetric(3),  rjb{i}   = dist_rjb4   (r0,rf,RW,RL,geom,ellipsoid);end
    if Rmetric(6),  rx{i}    = dist_rx4    (r0,rf,RW,RL,geom,ellipsoid);end
    if Rmetric(7),  ry0{i}   = dist_ry04   (r0,rf,RW,RL,geom);end
    if Rmetric(8),  zhyp{i}  = dist_zhyp4  (r0,rf,RW,RL,geom,ellipsoid);end
    if Rmetric(9),  ztor{i}  = dist_ztor4  (r0,rf,RW,RL,geom,ellipsoid);end
end

if Rmetric(1)==1
    M    = cell(size(rrup));
    rate  = cell(size(rrup));
    for i=1:nM
        nri =size(rrup{i},1);
        M{i}=MM(i)*ones(nri,1);
        rate{i}=1/nri*ones(nri,1)*rateM(i);
    end
else
    if Rmetric(2)==1
        M    = cell(size(rhyp));
        rate  = cell(size(rhyp));
        for i=1:nM
            nri =size(rhyp{i},1);
            M{i}=MM(i)*ones(nri,1);
            rate{i}=1/nri*ones(nri,1)*rateM(i);
        end
    elseif Rmetric(3)==1
        M    = cell(size(rjb));
        rate  = cell(size(rjb));
        for i=1:nM
            nri =size(rjb{i},1);
            M{i}=MM(i)*ones(nri,1);
            rate{i}=1/nri*ones(nri,1)*rateM(i);
        end
    end
end

M     = vertcat(M{:});
rate  = vertcat(rate{:});
if Rmetric(1),  rrup  = vertcat(rrup {:});end
if Rmetric(2),  rhyp  = vertcat(rhyp {:});end
if Rmetric(3),  rjb   = vertcat(rjb  {:});end
if Rmetric(6),  rx    = vertcat(rx   {:});end
if Rmetric(7),  ry0   = vertcat(ry0  {:});end
if Rmetric(8),  zhyp  = vertcat(zhyp {:});end
if Rmetric(9),  ztor  = vertcat(ztor {:});end

%% GMM PARAMETERS
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
    
    
    switch str_test
        case 'Youngs1997',              param = [M,rrup,zhyp,VS30,usp];
        case 'AtkinsonBoore2003',       param = [M,rrup,zhyp,VS30,usp];
        case 'Zhao2006',                param = [M,rrup,zhyp,VS30,usp];
        case 'Mcverry2006',             param = [M,rrup,zhyp,VS30,usp];
        case 'ContrerasBoroschek2012',  param = [M,rrup,zhyp,VS30,usp];
        case 'BCHydro2012',             param = [M,rrup,rhyp,zhyp,VS30,usp];
        case 'BCHydro2018',             param = [M,rrup,ztor,VS30,usp];
        case 'Kuehn2020',               param = [M,rrup,ztor,VS30,usp];
        case 'Parker2020',              param = [M,rrup,zhyp,VS30,usp];
        case 'Arteta2018',              param = [M,rhyp,VS30,usp];
        case 'Idini2016',               param = [M,rrup,rhyp,zhyp,VS30,usp];
        case 'MontalvaBastias2017',     param = [M,rrup,rhyp,zhyp,VS30,usp];
        case 'MontalvaBastias2017HQ',   param = [M,rrup,rhyp,zhyp,VS30,usp];
        case 'Montalva2018'
            [~,indf0] = intersect(hparam,'f0');    f0    = media(indf0);
            param = [M,rrup,rhyp,zhyp,VS30,f0,usp];
        case 'SiberRisk2019',           param = [M,rrup,rhyp,zhyp,VS30,usp];
        case 'Garcia2005',              param = [M,rrup,rhyp,zhyp,usp];
        case 'Jaimes2006',              param = [M,rrup,usp];
        case 'Jaimes2015',              param = [M,rrup,usp];
        case 'Jaimes2016',              param = [M,rrup,usp];
        case 'GarciaJaimes2017',        param = [M,rrup,usp];
        case 'GarciaJaimes2017HV',      param = [M,rrup,usp];
        case 'Bernal2014',              param = [M,rrup,zhyp,usp];
        case 'Sadigh1997',              param = [M,rrup,VS30,usp];
        case 'I2008',                   param = [M,rrup,VS30,usp];
        case 'CY2008',                  param = [M,rrup,rjb,rx,ztor,dip,VS30,usp];
        case 'BA2008',                  param = [M,rjb,VS30,usp];
        case 'CB2008',                  param = [M,rrup,rjb,ztor,dip,VS30,usp];
        case 'AS2008',                  param = [M,rrup,rjb,rx,ztor,dip,width,VS30,usp];
        case 'AS1997h',                 param = [M,rrup,VS30,usp];
        case 'I2014',                   param = [M,rrup,VS30,usp];
        case 'CY2014',                  param = [M,rrup,rjb,rx,ztor,dip,VS30,usp];
        case 'CB2014',                  param = [M,rrup,rjb,rx,zhyp,ztor,'unk',dip,width,VS30,usp];
        case 'BSSA2014',                param = [M,rjb,VS30,usp];
        case 'ASK2014',                 param = [M,rrup,rjb,rx,ry0,ztor,dip,width,VS30,usp];
        case 'AkkarBoomer2007',         param = [M,rjb,usp];
        case 'AkkarBoomer2010',         param = [M,rjb,usp];
        case 'Arroyo2010',              param = [M,rrup,usp];
        case 'Bindi2011',               param = [M,rjb,VS30,usp];
        case 'Kanno2006',               param = [M,rrup,zhyp,VS30,usp];
        case 'Cauzzi2015',              param = [M,rrup,rhyp,VS30,usp];
        case 'DW12',                    param = [M,rrup,usp];
        case 'FG15',                    param = [M,rrup,zhyp,VS30,usp];
        case 'TBA03',                   param = [M,rrup,usp];
        case 'BU17',                    param = [M,rrup,zhyp,usp];
        case 'CB10',                    param = [M,rrup,rjb,ztor,dip,VS30,usp];
        case 'CB11',                    param = [M,rrup,rjb,ztor,dip,VS30,usp];
        case 'CB19',                    param = [M,rrup,rjb,rx,zhyp,ztor,'unk',dip,width,VS30,usp];
        case 'KM06',                    param = [M,rrup,usp];
        case 'medianPCE_bchydro',       param = [M,rrup,VS30,pce];            

        case 'PCE_nga',                 param = [M,rrup,VS30,usp];
        case 'PCE_bchydro',             param = [M,rrup,VS30,usp];
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
                case 'Macedo2020', param = {M,rrup,VS30 ,gmm.cond,param{:}}; %#ok<CCAT>
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
