function data=layerdatainterp(p,layer,fd)

layer    = layer(end:-1:1);
baseline = layer{1};
if ischar(baseline)
    baseline = str2double(baseline);
end
Nsource  = length(layer);
Nsites   = size(p,1);
data     = ones(Nsites,1)*baseline;

if numel(layer)==1
    return
end

Lon = p(:,1);
Lat = p(:,2);

for i=2:Nsource
    DATA  = load(layer{i});
    fldname = fields(DATA);
    fldname = fldname{1};
    
    switch fldname
        case 'raster' % tiff file (e.g., from USGS)
            source = DATA.raster;
            Npoly = length(source);
            for j=1:Npoly
                r = source(j);
                IN  = inpolygon(Lat,Lon,r.box(:,2),r.box(:,1));
                if any(IN)
                    img       = double(imread(r.tif));
                    lon       = linspace(r.box(1,1),r.box(2,1),r.nlon);
                    lat       = linspace(r.box(2,2),r.box(1,2),r.nlat);
                    [LON,LAT] = meshgrid(lon,lat);
                    data(IN)  = interp2(LON,LAT,img,Lon(IN),Lat(IN));
                end
            end
        case 'microzone' %simple matfile
            source = DATA.microzone;
            Npoly = length(source);
            for j=1:Npoly
                r = source(j);
                IN  = inpolygon(Lon,Lat,r.Lon,r.Lat);
                if any(IN)
                    data(IN)  = source(j).(fd);
                end
            end
    end
end

