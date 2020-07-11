function[SOF]=rake2sof(rake)

% style of faulting from rake angle as defined in Aki and Richards (1980, p106)
if and(-30<=rake,rake<30) || and(150<=rake,rake<=180) || and(-180<=rake,rake<-150)
    SOF='strike-slip';
end

if and(30<=rake,rake<60) || and(120<=rake,rake<150)
    SOF='reverse-oblique';
end

if and(60<=rake,rake<120)
    SOF='reverse';
end

if and(-60<=rake,rake<-30) || and(-150<=rake,rake<-120)
    SOF='normal-oblique';
end

if and(-120<=rake,rake<-60)
    SOF='normal';
end

if rake==-999
    SOF='unspecified';
end
