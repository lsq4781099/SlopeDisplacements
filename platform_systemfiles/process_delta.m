function[obj]=process_delta(obj)

nObj = length(obj.source);
if nObj==0
    return
end

obj.mptr   = zeros(nObj,2);
obj.mdata  = zeros(nObj,2);
obj.meanMo = zeros(nObj,1);
for i=1:nObj
    M     = obj.num(i,3);
    Mchar = obj.num(i,3);
    mweight = 1;
    [mpdf,~,meanMo] = delta(M,Mchar);
    dPm = mweight.*mpdf/(mweight'*mpdf);
    obj.mptr(i,:) = [i,i];
    obj.mdata(i,:)=[M,dPm];
    obj.meanMo(i)=meanMo;
end

