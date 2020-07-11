function[obj]=createObj(tag)

switch tag
    % seismic source objects
    case 'model'
        obj  = struct('id',[],'source',[]);
        
    case 'source'
        obj  = struct(...
            'txt'      , [],...
            'adp'      , [],...
            'num'      , [],...
            'vptr'     , [],...
            'mptr'     , [],...
            'vert'     , [],...
            'center'   , [],...
            'xyzm'     , [],...
            'conn'     , [],...
            'aream'    , [],...
            'hypm'     , [],...
            'strike'   , [],...
            'normal'   , [],...
            'p'        , [],...
            'dsv'      , [],...
            'dip'      , [],...
            'mred'     , []);
        
    case 'sourceDSHA'
        obj  = struct(...
            'branch',[],...
            'txt'      , [],...
            'adp'      , [],...
            'num'      , [],...
            'vptr'     , [],...
            'mptr'     , [],...
            'vert'     , [],...
            'xyzm'     , [],...
            'conn'     , [],...
            'aream'    , [],...
            'hypm'     , [],...
            'strike'   , [],...
            'normal'   , [],...
            'p'        , [],...
            'dsv'      , [],...
            'dip'      , [],...
            'mred'     , []);
        
    case 'gmmlib'
        obj = struct('label',[],'txt',[],'type',[],'handle',[],'T',[],'usp',[],'Rmetric',[],'Residuals',[],'cond',[],'var',[]);
        
    case 'mscl'
        obj = struct('source',[],'adp',[],'num',[],'mptr',[],'mdata',[],'meanMo',[]);
        
    case 'site'
        obj.id        = cell(0,1);
        obj.p         = zeros(0,3);
        obj.t         = cell(0,2);
        obj.shape     = zeros(0,0);
        obj.param     = cell(0,1);
        obj.value     = zeros(0,0);
        
    case 'defaultlayers'
        obj.VS30{1}   = 760;      % Site Class
        obj.IdiniSiteClass{1}= 1; % Spectral Site class used in Idini.
        obj.To{1}     = 1;        % Fundamental Site Period
        obj.ky{1}     = 0.2;      % Yield Strength of Slops (used in PSDA2) 
        obj.covky{1}  = 0;        % Coef. of variation of Yield Strength of Slops (used in PSDA2)
        obj.Ts{1}     = 0.5;      % Mean slope period (used in PSDA2)
        obj.covTs{1}  = 0;        % Coef. of variation slope period (used in PSDA2)
        
    case 'opt'
        load pshatoolbox_RealValues opt
        obj=opt;
        
    case 'hazoptions'
        obj=struct('mod',1,'avg',[1 0 0 50 0],'sbh',[1 0 0 1],'dbt',[1 0 0 50],'map',[475 1],'pce',[0 1 50],'rnd',1);
        
    case 'Rbin'
        rmin  = 0;   rmax    = 360; dr    = 40;   obj  = [(rmin:dr:rmax-dr)',(rmin:dr:rmax-dr)'+dr];
        
    case 'Mbin'
        mmin  = 5;   mmax    = 9.6; dm    = 0.2;  obj  = [(mmin:dm:mmax-dm)',(mmin:dm:mmax-dm)'+dm];
        
    case 'Ebin'
        emin  = 0;   emax    = 2;   de    = 0.5;  obj  = [[-inf emin];(emin:de:emax-de)',(emin:de:emax-de)'+de;[emax,inf]];
        
    case 'returnperiods'
        obj=[144;250;475;949;1462;1950;2475;4975;10000];
        
    case 'LIBSsite'
        obj = struct('B',[],'L',[],'Df',[],'Q',[],'type',[],'wt',[],'LPC',[],'th1',[],'th2',[],'N1',[],'N2',[],'meth',[],'N160',[],'qc1N',[],'thick',[],'d2mat',[],'CPT',[]);
        
    case 'interfaceEQ'
        obj.Mech  = 1;
        obj.Mag   = 7;
        obj.Ztor  = 30;
        obj.Rrup  = 100;
        obj.Rx    = sqrt(obj.Rrup.^2-obj.Ztor.^2);
        obj.Rhyp  = obj.Rrup;
        obj.Zhyp  = obj.Ztor;
        
        Nsamples    = 40;
        obj.rpMag   = obj.Mag*ones(Nsamples,1);
        obj.rpZtor  = obj.Ztor*ones(Nsamples,1);
        obj.rpRrup  = logsp(obj.rpZtor(1),400,Nsamples)';
        obj.rpRx    = sqrt(max(obj.rpRrup.^2-obj.rpZtor.^2,0));
        obj.rpRhyp  = obj.rpRrup;
        obj.rpZhyp  = obj.rpZtor;
        
    case 'crustalEQ'
        obj.Mech   = 1;
        obj.Mag    = 7;
        obj.dip    = 90;
        obj.W      = 12;
        obj.Zbot   = 999;
        obj.Ztor   = 2;
        obj.Rx     = 100;
        obj        = getRrupRjb(obj);
        obj.Ry0    = 0;
        obj.region = 1;
        
        Nsamples   = 40;
        SC.Mag    = obj.Mech*ones(Nsamples,1);
        SC.dip    = obj.dip;
        SC.W      = obj.W;
        SC.Zbot   = obj.Zbot;
        SC.Ztor   = obj.Ztor*ones(Nsamples,1);
        SC.Rx     = [0;logsp(1,400,Nsamples-1)'];
        SC        = getRrupRjb(SC);
        SC.Ry0    = 0*ones(Nsamples,1);
        
        obj.rpMag  = SC.Mag;
        obj.rpZtor = SC.Ztor;
        obj.rpRx   = SC.Rx;
        obj.rpRrup = SC.Rrup;
        obj.rpRjb  = SC.Rjb;
        obj.rpRhyp = SC.Rhyp;
        obj.rpZhyp = SC.Zhyp;
        obj.rpRy0  = SC.Ry0;
        
    case 'siteGMM'
        obj.VS30   = 760;
        obj.f0     = 1;
        obj.Z10    = 0.048;
        obj.Z25    = 0.607;
        obj.Quali  = 'soil-soft';
        obj.Idini  = 2; %1=>sI,2=>sII,..., 6=>sVI
        obj.SGS    = 1; %1=>B,2=>C,3=>D
        
        
end