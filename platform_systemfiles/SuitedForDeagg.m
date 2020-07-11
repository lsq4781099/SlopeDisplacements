function dlist = SuitedForDeagg(handles)

if isempty(handles.sys.isREG)
    dlist=[];
else
    IJK       = handles.IJK;
    Nbranches = size(IJK,1);
    dlist     = 1:Nbranches;
    DD        = [];
    for i=1:Nbranches
        geomptr = handles.sys.branch(IJK(i,1),1);
        stype  = unique(handles.sys.mech{geomptr});
        [~,B]  = intersect([1;2;3],stype);
        usedM  = handles.T3(IJK(i,3),B+1);
        [~,B ] = intersect({handles.sys.SMLIB.id},usedM);
        I      = unique(vertcat(handles.sys.SMLIB(B).integrator));
        if any(ismember([3 4],I))
            DD = [DD;i]; %#ok<AGROW>
        end
    end
    dlist(DD)=[];
end