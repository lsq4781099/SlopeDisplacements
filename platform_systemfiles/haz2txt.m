function[]=haz2txt(FileName,PathName,WEIGHT,MRE,opt,h)

if isempty(MRE)
    return
end
id     = h.id;
p      = h.p;
im     = opt.im(:,1)';

Nsites = size(p,1);

FNAME = [PathName,FileName];
fileID = fopen(FNAME,'w');
NIM = length(opt.IM);

str   = IM2str(opt.IM);
Le = 0;
for i=1:length(str)
    Le = max(Le,1+length(str{i}));
end
% Le    = num2str(Le);
fprintf(fileID,' ------ Seismic Hazard Analysis ------\n');
fprintf(fileID,' Number of Sites   : %g\n',Nsites);
fprintf(fileID,' Hazard Points     : %g\n',length(im));
fprintf(fileID,' Max. Distance     : %g km\n',opt.MaxDistance);
fprintf(fileID,' ---------------------------------------\n');
fprintf(fileID,' \n');

text2 = repmat(' %12.8e',1,NIM);
text3 = '';
for i=1:NIM
    text3=[text3,sprintf('%-15s',str{i})]; %#ok<AGROW>
end

for site=1:Nsites
    LAMBDA  = nansum(MRE(site,:,:,:,:),4);
    LAMBDA  = permute(LAMBDA,[5,2,3,1,4]);
    lambda0 = prod(bsxfun(@power,LAMBDA,WEIGHT(:)),1);
    lambda0 = permute(lambda0,[3,2,1]);
    
    fprintf(fileID,'Site : %d \n',i);
    fprintf(fileID,'ID   : %s \n',id{site});
    if opt.ellipsoid.Code==0
        fprintf(fileID,'X(km): %6.4f\n',p(site,2));
        fprintf(fileID,'Y(km): %6.4f\n',p(site,1));
    else
        fprintf(fileID,'Lat  : %6.4f�\n',p(site,1));
        fprintf(fileID,'Lon  : %6.4f�\n',p(site,2));
    end
    fprintf(fileID,['im             ',text3,'\n']);
    fprintf(fileID,['%12.8e',text2,'\n'],[im;lambda0]);
    fprintf(fileID,' \n');
end
fclose(fileID);
if ispc,winopen(FNAME);end