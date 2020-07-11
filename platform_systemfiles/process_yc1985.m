function[obj]=process_yc1985(obj,options)

nObj  = length(obj.source);
if nObj==0
    return
end

switch options{1}
    case 'isampling'
        nM   = options{2};
        Mmin  = obj.num(:,4);
        Mchar = obj.num(:,5);
        Mmax  = Mchar+0.25;
        
        obj.mptr   = zeros(nObj,2);
        obj.mdata  = zeros(nM*nObj,2);
        obj.meanMo = zeros(nObj,1);
        alfa = 0.8;
        for i=1:nObj
            DM   = (Mmax(i)-Mmin(i))/sum(alfa.^(0:nM-1));
            M    = Mmin(i)+[0,DM*cumsum(alfa.^(0:nM-1))]';
            param   = obj.num(i,3:5);
            [~,mcdf,obj.meanMo(i)] = yc1985(M,param);
            M     = M(1:end-1)+diff(M)/2;
            dPm = diff(mcdf);
            ind = (1:nM)+nM*(i-1);
            obj.mptr(i,:)=ind([1 end]);
            obj.mdata(ind,:)=[M,dPm];
        end
        
        
    case 'gauss'
        nM0    = options{2};
        obj.mptr   = zeros(nObj,2);
        obj.mdata  = zeros(0,2);
        obj.meanMo = zeros(nObj,1);
        nM1   = max(round(nM0*4/5),5);
        nM2   = max(nM0-nM1,4);
        nM    = nM1+nM2;
        mshift=0;
        for i=1:nObj
            Mmin  = obj.num(i,4);
            Mchar = obj.num(i,5);
            mint1   = [Mmin,Mchar-0.25];
            mint2   = [Mchar-0.25,Mchar+0.25];
            [M1,mweight1] = gaussquad(nM1,mint1);
            [M2,mweight2] = gaussquad(nM2,mint2);
            M       = [M1;M2];
            mweight = [mweight1;mweight2];
            param   = obj.num(i,3:5);
            [mpdf,~,obj.meanMo(i)]  = yc1985(M,param);
            dPm = mweight.*mpdf/(mweight'*mpdf);
            obj.mptr(i,:)=[1 nM]+mshift;
            obj.mdata=[obj.mdata;[M,dPm]];
            mshift = size(obj.mdata,1);
        end
        
    case 'uniform'
        dM    = options{2};
        obj.mptr   = zeros(nObj,2);
        obj.mdata  = zeros(0,2);
        obj.meanMo = zeros(nObj,1);
        mshift=0;
        for i=1:nObj
            Mmin    = obj.num(i,4);
            Mmax    = obj.num(i,5)+0.25;
            M       = (Mmin+dM/2:dM:Mmax-dM/2)';
            nM      = length(M);
            mweight = dM*ones(size(M));
            param   = obj.num(i,3:5);
            [mpdf,~,obj.meanMo(i)]  = yc1985(M,param);
            dPm = mweight.*mpdf/(mweight'*mpdf);
            obj.mptr(i,:)=[1 nM]+mshift;
            obj.mdata=[obj.mdata;[M,dPm]];
            mshift = size(obj.mdata,1);
        end
end