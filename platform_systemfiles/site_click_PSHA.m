function site_click_PSHA(hObject,~,fig,p,HazMap,TriScatterd,entrymode)

ax1=findall(fig,'tag','ax1');
% ax2=findall(fig,'tag','ax2');

delete(findall(ax1,'Tag','satext'));

if entrymode==1
    gps    = [get(hObject,'XData')',get(hObject,'YData')'];
    point  = get(ax1,'CurrentPoint');
    point  = point(1,1:2);
    Ngps   = size(gps,1);
    select = point(ones(Ngps,1),:);
    disc   = sum((select-gps).^2,2);
    [~,B]  = min(disc);
    x=p(B,1);
    y=p(B,2);
    hazlevel = 1/HazMap.map(1);
    ch = findall(fig,'tag','lambda0');
    IMhaz = robustinterp(ch.YData,ch.XData,hazlevel,'loglog');
    if ~isnan(IMhaz)
        str = sprintf(' %3.2g',IMhaz);
        text(x,y,str,'parent',ax1,'Tag','satext','visible','on','verticalalignment','middle');
    end
end

if entrymode==2
    Tag = str2double(get(hObject,'Tag'));
    coordinates = get(ax1,'CurrentPoint');
    x   = coordinates(1,1);
    y   = coordinates(1,2);
    Fxy = TriScatterd{Tag}(x,y);
    str = sprintf(' %3.2g',Fxy);
    text(x,y,str,'parent',ax1,'Tag','satext','visible','on');
%     set(handles.shading,'ButtonDownFcn',{@site_click_PSHA,handles,2});
end


