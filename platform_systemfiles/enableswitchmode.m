function out=enableswitchmode(isMRE,isSHAKE)

if and(isMRE,isSHAKE)
    out='on';
else
    out='inactive';
end
