function plot_geometry_PSHA(ax,sys,opt)

vis='off';
if isempty(opt.ellipsoid.Code)
    xlabel(ax,'X (km)','fontsize',8,'fontname','arial','tag','xlabel')
    ylabel(ax,'Y (km)','fontsize',8,'fontname','arial','tag','ylabel')
else
    xlabel(ax,'Lon','fontsize',8,'fontname','arial','tag','xlabel')
    ylabel(ax,'Lat','fontsize',8,'fontname','arial','tag','ylabel')
end

meshcolor  = [0.6 0.6 1];
edgecolor   = [1 0.6 0.6];

Ngeom = size(sys.Nsrc,2);
for i=1:Ngeom
    n1 = length(sys.src1(i).txt);
    vert2=xyz2gps(sys.src2(i).xyzm , opt.ellipsoid); n2=size(vert2,1);
    vert3=xyz2gps(sys.src3(i).xyzm , opt.ellipsoid); n3=size(vert3,1);
    vert4=xyz2gps(sys.src4(i).xyzm , opt.ellipsoid); n4=size(vert4,1);
    vert5=xyz2gps(sys.src5(i).xyzm , opt.ellipsoid); n5=size(vert5,1);
    vert6=xyz2gps(sys.src6(i).xyzm , opt.ellipsoid); n6=size(vert6,1);
    vert7=xyz2gps(sys.src7(i).xyzm , opt.ellipsoid); n7=size(vert7,1);
    
    t2=[];
    if n7>0
        t  = sys.src7(i).conn;
        t2 = [t(:,[1 2 3]);t(:,[1 4 3]);t(:,[2 3 4]);t(:,[1 2 4])];
        t2 = unique(sort(t2,2),'rows');
    end
    
    % plot sourcemesh
    tag = sprintf('mesh%g',i); delete(findall(ax,'tag',tag));
    if n3>0,patch('parent',ax,'vertices',vert3(:,1:2),'faces',sys.src3(i).conn,'facecolor','none','edgecolor',meshcolor,'tag',tag,'visible',vis);end
    if n4>0,patch('parent',ax,'vertices',vert4(:,1:2),'faces',sys.src4(i).conn,'facecolor','none','edgecolor',meshcolor,'tag',tag,'visible',vis);end
    if n5>0,patch('parent',ax,'vertices',vert5(:,1:2),'faces',sys.src5(i).conn,'facecolor','none','edgecolor',meshcolor,'tag',tag,'visible',vis);end
    if n6>0,patch('parent',ax,'vertices',vert6(:,1:2),'faces',sys.src6(i).conn,'facecolor','none','edgecolor',meshcolor,'tag',tag,'visible',vis);end
    if n7>0,patch('parent',ax,'vertices',vert7(:,1:2),'faces',t2               ,'facecolor','none','edgecolor',meshcolor,'tag',tag,'visible',vis);end
    
    % plot sources
    tag = sprintf('edge%g',i); delete(findall(ax,'tag',tag));
    if n1>0,plot(ax,sys.src1(i).vert(:,1) ,sys.src1(i).vert(:,2) ,'ks','markersize',4,'markerfacecolor',edgecolor,'tag',tag,'visible',vis);end
    if n2>0,plot(ax,sys.src2(i).vert(:,1) ,sys.src2(i).vert(:,2) ,'color',edgecolor,'tag',tag,'linewidth',1,'visible',vis);end
    if n3>0,plot(ax,sys.src3(i).vert(:,1) ,sys.src3(i).vert(:,2) ,'color',edgecolor,'tag',tag,'visible',vis);end
    if n4>0,plot(ax,sys.src4(i).vert(:,1) ,sys.src4(i).vert(:,2) ,'color',edgecolor,'tag',tag,'visible',vis);end
    if n5>0,plot(ax,sys.src5(i).vert(:,1) ,sys.src5(i).vert(:,2) ,'color',edgecolor,'tag',tag,'visible',vis);end
    if n6>0,plot(ax,sys.src6(i).vert(:,1) ,sys.src6(i).vert(:,2) ,'color',edgecolor,'tag',tag,'visible',vis);end
    if n7>0,plot(ax,sys.src7(i).vert(:,1) ,sys.src7(i).vert(:,2) ,'color',edgecolor,'tag',tag,'visible',vis);end
    
    tag = sprintf('mesh%g',i);
    if n2>0, plot(ax,vert2(:,1),vert2(:,2),'k.','tag',tag,'visible',vis); end
    
    % source tag
    tag=sprintf('sourcetag%g',i);
    txt = [
        sys.src2(i).txt;
        sys.src3(i).txt;
        sys.src4(i).txt;
        sys.src5(i).txt;
        sys.src6(i).txt;
        sys.src7(i).txt];
    cm  = [
        sys.src2(i).center;
        sys.src3(i).center;
        sys.src4(i).center;
        sys.src5(i).center;
        sys.src6(i).center
        sys.src7(i).center];
    
    text(ax,cm(:,1),cm(:,2),strrep(txt,'_','-'),'fontsize',8,'fontweight','bold','visible',vis,'tag',tag,'horizontalalignment','center','verticalalignment','bottom');
end

