function[pages]=setDSHApages(Nscen)

nMax   = min(Nscen,20);
pages  = unique([1:nMax:Nscen,Nscen])';
if numel(pages)==1
    pages = [1 1];
else
    pages  = [pages(1:end-1),pages(2:end)-1];pages(end)=Nscen;
end