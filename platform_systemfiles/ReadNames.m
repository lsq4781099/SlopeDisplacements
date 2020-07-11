function [ FList ] = ReadNames(DataFolder)

DirContents=dir(DataFolder);
FList=[];
if(strcmpi(computer,'PCWIN') || strcmpi(computer,'PCWIN64'))
    NameSeperator='\';
elseif(strcmpi(computer,'GLNX86') || strcmpi(computer,'GLNXA86'))
    NameSeperator='/';
end

for i=1:numel(DirContents)
    if(~(strcmpi(DirContents(i).name,'.') || strcmpi(DirContents(i).name,'..')))
        if(~DirContents(i).isdir)
            [~,~,extension]=fileparts(DirContents(i).name);
            if any(strcmp(extension,{'.fig','.txt','.m','.mat','.png','.jpg'}))
                FList=cat(1,FList,{[DataFolder,NameSeperator,DirContents(i).name]});
            end
        else
            getlist=ReadNames([DataFolder,NameSeperator,DirContents(i).name]);
            FList=cat(1,FList,getlist);
        end
    end
end
end
