function[indx,hc]=compute_clusters(opt,h)

indx = [];
hc   = createObj('site');

switch opt.Clusters{1}
    case 'off'
        
    case 'on'
        Nsites = length(h.id);
        if Nsites==0
            return
        end
        
        Nc = min(opt.Clusters{2}(1),Nsites);
        if Nsites==Nc
            indx = 1:Nsites;
            hc   = h;
        else
            Nr = opt.Clusters{2}(2);
            stream  = RandStream('mlfg6331_64');  % Random number stream
            options = statset('UseParallel',1,'UseSubstreams',1,'Streams',stream);
            
            % Lat,Lon,Elev clusters
            [indx,Y] = kmeans([h.p,h.value],Nc,'Options',options,'MaxIter',10000,'Display','final','Replicates',Nr);
            hc.id    = compose('C%i',1:Nc)';
            hc.p     = Y(:,1:3);
            hc.param = h.param;
            hc.value = Y(:,4:end);
        end
        
end
