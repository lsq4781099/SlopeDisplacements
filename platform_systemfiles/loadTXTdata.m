function[data]=loadTXTdata(data,tgtxt)

ellip = referenceEllipsoid('WGS84','km');

ind   = strfind(data,tgtxt);
geom = [];
for i=1:size(data,1)
    if ~isempty(ind{i}), geom=[geom;i]; end %#ok<*AGROW>
end

if isempty(geom)
    return
end

list0 = (1:size(data,1))';
list0(geom)=-1;

%% add new geom
Ngeom = length(geom);
t(1:Ngeom)        = cell(Ngeom,1);
newlist1(1:Ngeom) = cell(Ngeom,1);
for jj=1:Ngeom
    str = regexp(data{geom(jj)},'\ ','split');
    filename = str{2};
    
    % preprocess input string
    do_search = false;
    if ~strcmp(str{4},'all')
        R = str2double(str(4));
        center = str2double(str(5:6));
        do_search=true;
    end
    
    % variables to overwrite
    [~,B]=intersect(str,'overwrite');
    if length(str)>B
        overwrite = struct(str{B+1:end});
    else
        overwrite = [];
    end
    
    % reads file
    fid     = fopen(filename);
    line    = regexp(fgetl(fid),'\ ','split');
    if numel(line)==2 && do_search
        fname = line{2};
        load(fname,'loc')
        xyz     = gps2xyz(loc,ellip);
        xyz0    = gps2xyz(center,ellip);
        IND     = sum((xyz-xyz0).^2,2)<=R^2;
        Npoints = sum(IND);
        ptrs    = find(IND);
        skip    = diff([0;ptrs(:)])-1;
        t{jj}   = cell(Npoints,1);
        for i=1:Npoints
            line1    = textscan(fid,'%s',1,'headerLines',skip(i),'delimiter','\n');
            t{jj}{i} = line1{1}{1};
        end
        t{jj}=regexprep(t{jj},' +',' ');
    else
        t{jj} = textscan(fid,'%s','delimiter','\n');
        t{jj}=t{jj}{1};
        Npoints=length(t{jj});
    end
    
    % text overwrite
    if~isempty(overwrite)
        flds   = fields(overwrite);
        Nfield = length(flds);
        for kk=1:Nfield
            t{jj}=strrep(t{jj},sprintf('$%s',flds{kk}),overwrite.(flds{kk}));
        end
    end
    
    newlist1{jj} = geom(jj)+(1:Npoints)'/(Npoints+1);
    fclose(fid);
    
end

data(list0<0)=[];
list0(list0<0)=[];
for jj=1:Ngeom
    data=[data;strtrim(regexprep(t{jj},' +',' '))];
    t{jj}=[];
    list0=[list0;newlist1{jj}];
end
[~,sindex]=sort(list0);
data = data(sindex);

