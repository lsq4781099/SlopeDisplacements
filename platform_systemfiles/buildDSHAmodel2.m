function[source]=buildDSHAmodel2(sys,opt,h,bptr,sptr)

if nargin==3
    branch = sys.branch(1,:);
else
    branch = sys.branch(bptr,:);
end

geomptr   = branch(1);
gmpeptr   = branch(2);
msclptr   = branch(3);
epistemic = branch(4);

src1      = sys.src1(geomptr);
src2      = sys.src2(geomptr);
src3      = sys.src3(geomptr);
src4      = sys.src4(geomptr);
src5      = sys.src5(geomptr);
src6      = sys.src6(geomptr);
src7      = sys.src7(geomptr);
gmmptr    = sys.gmmptr(gmpeptr,:);
gmmlib    = sys.gmmlib(gmmptr);

mrr1     = sys.mrr1(msclptr);
mrr2     = sys.mrr2(msclptr);
mrr3     = sys.mrr3(msclptr);
mrr4     = sys.mrr4(msclptr);
mrr5     = sys.mrr5(msclptr);

numG  = sys.numG{geomptr};
numM  = sys.numM{msclptr};
cgm   = sys.cgm{geomptr,msclptr};
Nsrc   = size(cgm,1);


if nargin==3
    sptr = 1:size(numG,1);
else
    sptr = sptr(:)';
end

% assembly of area sources
source(1:Nsrc,1)=createObj('sourceDSHA');
for j=sptr
    source(j,1).branch=branch;
    spt =numG(j,2);
    switch numG(j,1)
        case 1
            source(j).obj     = 1;
            source(j).sourcelabel = src1.txt{spt};
            source(j).numgeom = src1.num(spt,:);
            source(j).hypm    = src1.vert(spt,:);
            source(j).aream   = 1;
            source(j).normal  = src1.normal(spt,:);
            source(j).gmm     = gmmlib(src1.num(spt,5));
            
        case 2
            mptr = src2.mptr(spt,:);
            ind  = mptr(1):mptr(2);
            source(j).obj     = 2;
            source(j).sourcelabel = src2.txt{spt};
            source(j).numgeom = src2.num(spt,:);
            source(j).hypm    = src2.hypm(ind,:);
            source(j).aream   = src2.aream(ind);
            source(j).normal  = src2.normal(ind,:);
            source(j).gmm     = gmmlib(src2.num(spt,5));
            
        case 3
            vptr = src3.vptr(spt,:);
            mptr = src3.mptr(spt,:);
            ind  = mptr(1):mptr(2);
            source(j).obj     = 3;
            source(j).sourcelabel = src3.txt{spt};
            source(j).numgeom = src3.num(spt,:);
            source(j).hypm    = src3.hypm(ind,:);
            source(j).aream   = src3.aream(ind);
            source(j).strike  = src3.strike(ind,:);
            source(j).normal  = src3.normal(ind,:);
            source(j).gmm     = gmmlib(src3.num(spt,5));
            source(j).vert    = src3.vert(vptr(1):vptr(2),:);
            
        case 4
            mptr = src4.mptr(spt,:);
            ind  = mptr(1):mptr(2);
            source(j).obj     = 4;
            source(j).sourcelabel = src4.txt{spt};
            source(j).numgeom = src4.num(spt,:);
            source(j).hypm    = src4.hypm(ind,:);
            source(j).aream   = src4.aream(ind);
            source(j).normal  = src4.normal(ind,:);
            source(j).gmm     = gmmlib(src4.num(spt,5));
            
            
        case 5
            mptr = src5.mptr(spt,:);
            ind  = mptr(1):mptr(2);
            source(j).obj     = 5;
            source(j).sourcelabel = src5.txt{spt};
            source(j).numgeom = src5.num(spt,:);
            source(j).hypm    = src5.hypm(ind,:);
            source(j).aream   = src5.aream(ind);
            source(j).p       = src5.p(:,:,ind);
            source(j).dsv     = src5.dsv(:,:,ind);
            source(j).gmm     = gmmlib(src5.num(spt,5));
            source(j).normal  = -src5.dsv(:,3)';
            source(j).dip     = source(j).numgeom(7);
            
            
        case 6
            vptr = src6.vptr(spt,:);
            mptr = src6.mptr(spt,:);
            
            ind  = mptr(1):mptr(2);
            reduction         = src6.num(spt,6);
            
            source(j).obj     = 6;
            source(j).sourcelabel = src6.txt{spt};
            source(j).numgeom = src6.num(spt,:);
            source(j).xyzm    = src6.xyzm;
            source(j).conn    = src6.conn(ind,:);
            source(j).hypm    = src6.hypm(ind,:);
            source(j).aream   = src6.aream(ind);
            source(j).normal  = src6.normal(ind,:);
            source(j).gmm     = gmmlib(src6.num(spt,5));
            source(j).vert    = src6.vert(vptr(1):vptr(2),:);
            
            % random permutation
            n = length(ind);
            k = round(n*reduction/100);
            IND = sort(randperm(n,k));
            source(j).aream  = source(j).aream(IND,:);
            source(j).hypm   = source(j).hypm(IND,:);
            source(j).normal = source(j).normal(IND,:);
            
        case 7
            mptr = src7.mptr(spt,:);
            ind  = mptr(1):mptr(2);
            Nind = length(ind);
            source(j).obj     = 7;
            source(j).numgeom = src7.num(spt,:);
            source(j).hypm    = src7.hypm(ind,:);
            source(j).aream   = src7.aream(ind);
            source(j).normal  = repmat(src7.normal,Nind,1);
            source(j).gmm     = gmmlib(src7.num(spt,5));
    end
    
    % updats usp if mechanism is read from source definition
    switch source(j).gmm.type
        case {'regular','pce'}
            [A,B]=ismember(source(j).gmm.txt,{'auto'});
            if any(A)
                rake = source(j).numgeom(2);
                vv   = find(B)-3;
                source(j).gmm.usp{vv} = rake2sof(rake);
            end
        case 'cond'
            [A,B]=ismember(source(j).gmm.txt,{'auto'});
            if any(A)
                rake = source(j).numgeom(2);
                vv   = find(B)-6;
                source(j).gmm.usp{vv} = rake2sof(rake);
            end
        case 'frn'
            for uu=1:length(source(j).gmm.usp)
                [A,B]=ismember(source(j).gmm.usp{uu}.txt,{'auto'});
                if any(A)
                    rake = source(j).numgeom(2);
                    vv   = find(B)-3;
                    source(j).gmm.usp{uu}.usp{vv} = rake2sof(rake);
                end
            end
    end
    
    % magnitude recurrence relation
    mp    = cgm(j);
    mtype = numM(mp,1);
    mptr  = numM(mp,2);
    switch mtype
        case 1
            source(j).nummag = mrr1.num(mptr,:);
            source(j).mscl   = mrr1.mdata(mptr,:);
            source(j).meanMo = mrr1.meanMo(mptr,:);
        case 2
            source(j).nummag = mrr2.num(mptr,:);
            ind = mrr2.mptr(mptr,1):mrr2.mptr(mptr,2);
            source(j).mscl   = mrr2.mdata(ind,:);
            source(j).meanMo = mrr2.meanMo(mptr,:);
        case 3
            source(j).nummag = mrr3.num(mptr,:);
            ind = mrr3.mptr(mptr,1):mrr3.mptr(mptr,2);
            source(j).mscl   = mrr3.mdata(ind,:);
            source(j).meanMo = mrr3.meanMo(mptr,:);
        case 4
            source(j).nummag = mrr4.num(mptr,:);
            ind = mrr4.mptr(mptr,1):mrr4.mptr(mptr,2);
            source(j).mscl   = mrr4.mdata(ind,:);
            source(j).meanMo = mrr4.meanMo(mptr,:);
        case 5
            source(j).nummag = mrr5.num(mptr,:);
            ind = mrr5.mptr(mptr,1):mrr5.mptr(mptr,2);
            source(j).mscl   = mrr5.mdata(ind,:);
            source(j).meanMo = mrr5.meanMo(mptr,:);
    end
end

for j=sptr
    if source(j).nummag(1)==1
        source(j).NMmin=source(j).nummag(2);
    elseif source(j).nummag(1)==2
        A          = sum(source(j).aream)* 1e10;  % cm2
        S          = source(j).nummag(2)*0.100;   % cm/year
        meanMo     = source(j).meanMo;
        source(j).NMmin = opt.ShearModulus*A*S/meanMo;
    end
end

source = source(sptr);

xyz = gps2xyz(h.p,opt.ellipsoid);
Nsite  = size(xyz,1);
ind    = zeros(Nsite,length(source));
for i=1:Nsite
    ind(i,:)=selectsource(opt.MaxDistance,xyz(i,:),source);
end

if ~isempty(source(sum(ind,1)>0))
    source=source(sum(ind,1)>0);
end


