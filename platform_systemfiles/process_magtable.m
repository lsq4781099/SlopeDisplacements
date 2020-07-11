function[obj]=process_magtable(obj)

nObj  = length(obj.source);
if nObj==0
    return
end
obj.mptr   = zeros(nObj,2);
obj.mdata  = zeros(0,2);
obj.meanMo = zeros(nObj,1);
mshift=0;
for i=1:nObj
    minMag   = obj.num{i}(3);
    binWidth = obj.num{i}(4);
    lambdaM  = obj.num{i}(5:end)';
    Nm       = length(lambdaM);
    Mmin     = minMag;
    Mmax     = (Mmin+binWidth*(Nm-1));
    M        = (Mmin:binWidth:Mmax)';
    dPm      = lambdaM/sum(lambdaM);
    nM       = length(M);
    obj.mptr(i,:)=[1 nM]+mshift;
    obj.mdata = [obj.mdata;[M,dPm]];
    obj.meanMo(i)= nan;
    mshift = size(obj.mdata,1);
end