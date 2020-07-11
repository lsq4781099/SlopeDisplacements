function[obj]=process_trexpub(obj,options)

nObj  = length(obj.source);
if nObj==0
    return
end

WORK IN PROGRESS

switch options{1}
    case 'isampling'
        nM   = options{2};
        Mmin = obj.num(:,4);
        Mmax = obj.num(:,5);
        obj.mptr   = zeros(nObj,2);
        obj.mdata  = zeros(nM*nObj,2);
        obj.meanMo = zeros(nObj,1);
        alfa = 0.8;
        for i=1:nObj
            DM   = (Mmax(i)-Mmin(i))/sum(alfa.^(0:nM-1));
            M    = Mmin(i)+[0,DM*cumsum(alfa.^(0:nM-1))]';
            param = obj.num(i,3:6);
            [~,mcdf,obj.meanMo(i)] = truncexp(M,param);
            M     = M(1:end-1)+diff(M)/2;
            dPm = diff(mcdf);
            ind = (1:nM)+nM*(i-1);
            obj.mptr(i,:)=ind([1 end]);
            obj.mdata(ind,:)=[M,dPm];
        end
        
    case 'gauss'
        nM   = options{2};
        Mmin = obj.num(:,4);
        Mmax = obj.num(:,5);
        obj.mptr   = zeros(nObj,2);
        obj.mdata  = zeros(nM*nObj,2);
        obj.meanMo = zeros(nObj,1);
        for i=1:nObj
            minterval   = [Mmin(i),Mmax(i)];
            [M,mweight] = gaussquad(nM,minterval);
            param = obj.num(i,3:6);
            [mpdf,~,obj.meanMo(i)] = truncexp(M,param);
            dPm = mweight.*mpdf/(mweight'*mpdf);
            ind = (1:nM)+nM*(i-1);
            obj.mptr(i,:)=ind([1 end]);
            obj.mdata(ind,:)=[M,dPm];
        end
    case 'uniform'
        dM   = options{2};
        Mmin = obj.num(:,4);
        Mmax = obj.num(:,5);
        obj.mptr   = zeros(nObj,2);
        obj.mdata  = zeros(0,2);
        obj.meanMo = zeros(nObj,1);
        mshift     = 0;
        for i=1:nObj
            M  = (Mmin(i)+dM/2:dM:Mmax(i)-dM/2)';
            nM = length(M);
            mweight = dM*ones(size(M));
            param = obj.num(i,3:6);
            [mpdf,~,obj.meanMo(i)] = truncexp(M,param);
            dPm = mweight.*mpdf/(mweight'*mpdf);
            obj.mptr(i,:)=[1 nM]+mshift;
            obj.mdata    = [obj.mdata;[M,dPm]];
            mshift = size(obj.mdata,1);
        end        
        
        
end
