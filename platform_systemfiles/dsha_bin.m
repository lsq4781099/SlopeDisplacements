
switch etype
    case 'IS'
        [rate,Y] =dsha_is(handles);
        fid = fopen(FILE,'w');
        data = [rate,Y];
        data = [data;ones(1,size(data,2))*(-999)];
        fwrite(fid,data(:),'double');
        fclose(fid);
        fprintf('done\n')
        
    case 'KM'
        if isempty(handles.krate) && isempty(handles.kY)
            [rate,Y] =dsha_kmeans(handles,handles.optkm);
        else
            rate = handles.krate;
            Y    = handles.kY;
        end
        fid = fopen(FILE,'w');
        data = [rate,Y];
        data = [data;ones(1,size(data,2))*(-999)];
        fwrite(fid,data,'double');
        fclose(fid);
        fprintf('done\n')
end