function tlist=UHSperiods(handles)

Tmax = [];
gptr = unique(unique(handles.sys.gmmptr(:)))';
for i=gptr
    Tmax = [Tmax;max(handles.sys.gmmlib(i).T)]; %#ok<AGROW>
end
Tmax  = min(Tmax);
tlist = unique([logsp(0.01,10,30)';Tmax]);
tlist(tlist>Tmax)=[];
