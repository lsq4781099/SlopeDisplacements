function u=getGRDfield(fig,sys,opt,data,h,uhs,modo)

% mod = 0, gets fields from PSHA
% mod = 1, gets fields from user defined UHS
% GRD analysis.
% This code performs a super-fast search of IM associated to the
% specified return period.
branch       = sys.branch(:,1:3);
ShearModulus = opt.ShearModulus;

switch modo
    case 'PSHA'
        val       = get(findall(fig,'tag','select_IM'),'value');
        model_ptr = get(findall(fig,'tag','IM_menu'),'value');
        site_ptr  = get(findall(fig,'tag','site_menu'),'value');
        IM        = data.T(val-1);
        h.id      = h.id(site_ptr);
        h.p       = h.p(site_ptr,:);
        h.value   = h.value(site_ptr,:);
        haz       = 1/str2double(get(findall(fig,'tag','retperiod'),'string'));
        RSR       = data.MapRSR(:,:,val-1);
        
        sources   = buildmodelin(sys,branch(model_ptr,:),ShearModulus);
        Ns        = length(sources);
        
        % round 1: seach of PGA
        im        = logsp(0.01,3,5)';
        MRE       = runhazard1(im,IM,h,opt,sources,Ns,1);
        MRE       = nansum(MRE,4)';
        IM0       = exp(interp1(log(MRE),log(im),log(haz),'spline'));
        
        % round 2
        immin     = 0.98*IM0;
        immax     = 1.02*IM0;
        im        = logsp(immin,immax,5)';
        MRE       = runhazard1(im,IM,h,opt,sources,Ns,1);
        MRE       = nansum(MRE,4)';
        IM0       = exp(interp1(log(MRE),log(im),log(haz),'spline'));
        u         = RSR*IM0;
    case 'Event'
        val       = get(findall(fig,'tag','select_IM'),'value');
        To        = data.T(val-1);
        x         = uhs(:,1);
        y         = uhs(:,2);
        IM0       = exp(interp1(log(x),log(y),log(To),'spline'));
        RSR       = data.MapRSR(:,:,val-1);
        u         = RSR*IM0;
end
end
