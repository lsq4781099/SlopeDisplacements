function[Data]=sSCEN_list(source)

Nsource = length(source);
Data=zeros(0,9);
for i=1:Nsource
    if any(source(i).obj==[1 2 3 4 7])
        nR = size(source(i).hypm,1);
        nM = size(source(i).mscl,1);
        [JJ,II]=meshgrid(1:nM,1:nR);
        II=II(:);
        JJ=JJ(:);
        
        id      = i*ones(nR*nM,1);
        M       = source(i).mscl(JJ,1);
        hypm    = source(i).hypm(II,:);
        rateM   = source(i).mscl(JJ,2);
        rateR   = source(i).aream(II,:)/sum(source(i).aream);
        rate    = source(i).NMmin*rateM.*rateR;
        normal  = source(i).normal(II,:);
        newData = [id,M,hypm,normal,rate];
        Data    = [Data;newData]; %#ok<AGROW>
    elseif source(i).obj==6
        nR = size(source(i).hypm,1);
        nM = size(source(i).mscl,1);
        [JJ,II]=meshgrid(1:nM,1:nR);
        II=II(:);
        JJ=JJ(:);
        
        id      = i*ones(nR*nM,1);
        M       = source(i).mscl(JJ,1);
        hypm    = source(i).hypm(II,:);
        rateM   = source(i).mscl(JJ,2);
        rateR   = source(i).aream(II,:)/sum(source(i).aream);
        rate    = source(i).NMmin*rateM.*rateR;
        normal  = source(i).normal(II,:);
        newData = [id,M,hypm,normal,rate];
        Data    = [Data;newData]; %#ok<AGROW>        
    elseif source(i).obj==5
        
        Area    = source(i).aream;
        MM      = source(i).mscl(:,1);
        rateM   = source(i).mscl(:,2);
        RA      = source(i).numgeom(3:4);
        L       = sum(abs(source(i).p(:,1)))/2;%source.numgeom(8);

        if source(i).numgeom(9)==-999
            width = Area/L;
        else
            width   = source(i).numgeom(9);
        end
        
        aratio  = source(i).numgeom(10);
        spacing = source(i).numgeom(11);
        
        % RUPTURE AREA AND SCENARIOS
        if RA(2)==0
            rupArea   = rupRelation(MM,0,RA(1));
        else
            NA = 25;
            x  = linspace(-2,2,NA)';
            dx = x(2)-x(1);
            rateA = normcdf(x+dx/2,0,1)-normcdf(x-dx/2,0,1);
            rateA = rateA/sum(rateA);
            rupArea = zeros(length(MM),NA);
            for ii=1:NA
                rupArea(:,ii) = rupRelation(MM,x(ii),RA(1));
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
        p     = source(i).p;
        pmean = source(i).hypm;
        rot   = source(i).dsv;
        rf    = cell(1,nM);
        
        for jj=1:nM
            RL    = rupLength(jj);
            RW    = rupWidth (jj);
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
            locxy = [locx(:),locy(:),zeros(nR,1)];
            rf{jj}      = pmean+locxy*rot';
            
        end
        
        M    = cell(size(rf));
        rate  = cell(size(rf));
        for j=1:nM
            nri =size(rf{j},1);
            M{j}=MM(j)*ones(nri,1);
            rate{j}=1/nri*ones(nri,1)*rateM(j);
        end
        M       = vertcat(M{:});
        hypm    = vertcat(rf{:});
        nScen   = length(M);
        normal  = repmat(-rot(:,3)',nScen,1);
        rate    = source(i).NMmin*vertcat(rate{:});
        id      = i*ones(nR*nM,1);
        newData = [id,M,hypm,normal,rate];
        Data    = [Data;newData]; %#ok<AGROW>
    end
end
