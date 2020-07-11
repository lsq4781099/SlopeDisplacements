function plot_hazard_PSHA(fig,haz,opt,isPCE,MRE,MREPCE,validation,weights,labelG,mech,idx,holdcolorbar)
            
ax2          = findall(fig,'tag','ax2');
addLeg       = findall(fig,'tag','addLeg');
ExportHazard = findall(fig,'tag','ExportHazard');

if isempty(MRE)
    return
end
ch=findall(ax2,'Type','line'); delete(ch);

if nargin ==11
   holdcolorbar=false;
end

if holdcolorbar==0
    ch=findall(fig,'tag','Colorbar'); delete(ch);
end
ch=findall(fig,'type','legend');delete(ch);


switch haz.mod
    case 1, str = plotHazMode1(fig,haz,opt,isPCE,MRE,MREPCE,validation,weights,idx); % plots averaged seismic hazard
    case 2, str = plotHazMode2(fig,haz,opt,MRE,labelG,mech,idx);                     % plots single branch hazard for GMMs of type REGULAR
    case 3, str = plotHazMode3(fig,haz,opt,MREPCE,idx);                              % plots single branch hazard for GMMs of type PCE
end

Leg=legend(ax2,strrep(str,'_',' '));
switch addLeg.Value
    case 0,Leg.Visible='off';
    case 1,Leg.Visible='on';
end
Leg.FontSize=8;
Leg.EdgeColor=[1 1 1];
Leg.Location='SouthWest';
Leg.Tag='hazardlegend';
set(ExportHazard,'enable','on');
