function[obj]=process_catalog(obj,source,options)

nObj  = length(obj.source);
if nObj==0
    return
end


[~,~,ext]= fileparts(obj.adp{1});
switch ext
    case '.mat', load(obj.adp{1});
    case '.xlsx', cat = importdata(obj.adp{1});
end

param(1:nObj) = struct('Mmin',[],'Mmax',[],'NMmin',[],'bvalue',[],'sigma_a',[],'sigma_b',[]);
for i=1:nObj
    vptr  = source.vptr(i,1):source.vptr(i,2);
    FX    = source.vert(vptr,1)';
    FY    = source.vert(vptr,2)';
    param(i) = runWeichert(obj.num(i,3:4),FX,FY,cat);
end

nM   = options{2};
Mmin = vertcat(param.Mmin);
Mmax = vertcat(param.Mmax);
bvalue = vertcat(param.bvalue);
NMmin = vertcat(param.NMmin);
obj.num(:,2)=NMmin;
param = [bvalue,Mmin,Mmax];
switch options{1}
    case 'isampling'
        obj.mptr   = zeros(nObj,2);
        obj.mdata  = zeros(nM*nObj,2);
        obj.meanMo = zeros(nObj,1);
        alfa = 0.8;
        for i=1:nObj
            DM   = (Mmax(i)-Mmin(i))/sum(alfa.^(0:nM-1));
            M    = Mmin(i)+[0,DM*cumsum(alfa.^(0:nM-1))]';
            [~,mcdf,obj.meanMo(i)] = truncexp(M,param(i,:));
            M     = M(1:end-1)+diff(M)/2;
            dPm = diff(mcdf);
            ind = (1:nM)+nM*(i-1);
            obj.mptr(i,:)=ind([1 end]);
            obj.mdata(ind,:)=[M,dPm];
        end
        
    case 'gauss'
        obj.mptr   = zeros(nObj,2);
        obj.mdata  = zeros(nM*nObj,2);
        obj.meanMo = zeros(nObj,1);
        for i=1:nObj
            minterval   = [Mmin(i),Mmax(i)];
            [M,mweight] = gaussquad(nM,minterval);
            [mpdf,~,obj.meanMo(i)] = truncexp(M,param(i,:));
            dPm = mweight.*mpdf/(mweight'*mpdf);
            ind = (1:nM)+nM*(i-1);
            obj.mptr(i,:)=ind([1 end]);
            obj.mdata(ind,:)=[M,dPm];
        end
    case 'uniform'
        obj.mptr   = zeros(nObj,2);
        obj.mdata  = zeros(0,2);
        obj.meanMo = zeros(nObj,1);
        mshift     = 0;
        for i=1:nObj
            M  = (Mmin(i)+dM/2:dM:Mmax(i)-dM/2)';
            nM = length(M);
            mweight = dM*ones(size(M));
            [mpdf,~,obj.meanMo(i)] = truncexp(M,param(i,:));
            dPm = mweight.*mpdf/(mweight'*mpdf);
            obj.mptr(i,:)=[1 nM]+mshift;
            obj.mdata    = [obj.mdata;[M,dPm]];
            mshift = size(obj.mdata,1);
        end        
        
        
end
