function plotLayers(handles)

val = handles.layerpop.Value;
fld = handles.layerpop.String{val};
layer  = handles.layer.(fld);
delete(findall(handles.ax1,'tag','layer'))
if length(layer)==1
    return
end

VIS = handles.DisplayMicrozones.Checked;


for jj = 1:length(layer)-1
    fname = layer{jj};
    DATA  = load(fname);
    
    % microzonation data
    if isfield(DATA,'microzone')
        
        coldata    = vertcat(DATA.microzone.(fld));
        ccluster   = kmeans(coldata,min(length(unique(coldata)),length(coldata)));
        Ncol       = max(ccluster);
        lineStyles = jet(Ncol);
        for i=1:length(DATA.microzone)
            LAT=DATA.microzone(i).Lat;LAT(isnan(LAT))=[];
            LON=DATA.microzone(i).Lon;LON(isnan(LON))=[];
            patch(handles.ax1,LON,LAT,lineStyles(ccluster(i),:),'tag','layer','visible',VIS)
        end
    end
    
    % tiff data
    if isfield(DATA,'raster')
        for i=1:length(DATA.raster)
            LON = DATA.raster(i).box([1 2 2 1],1);
            LAT = DATA.raster(i).box([1 1 2 2],2);
            patch(handles.ax1,LON,LAT,[0.92941       0.6902      0.12941],'facealpha',0.2,'tag','layer','visible',VIS,'edgecolor','none')
        end
    end
    
end

