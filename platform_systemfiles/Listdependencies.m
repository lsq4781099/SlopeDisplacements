clc

% fun = 'SeismicHazard.m';
fun = 'GMselect.m';
% fun = 'PSDA2.m';
% fun = 'pshatoolbox_methods.m';
list  = matlab.codetools.requiredFilesAndProducts(fun)';
for i=1:size(list,1)
    copyfile(list{i},'C:\Users\Gabriel\Desktop\GMSelect')
end

%% Pcode
D=dir('*.m'); for i=1:size(D,1); pcode(D(i).name); end



%%
clearvars
clc
fun = 'SeismicHazard.m';        list1  = matlab.codetools.requiredFilesAndProducts(fun)';1
fun = 'GMselect.m';             list2  = matlab.codetools.requiredFilesAndProducts(fun)';2
fun = 'PSETTLE.m';              list3  = matlab.codetools.requiredFilesAndProducts(fun)';3
fun = 'pshatoolbox_methods.m';  list4  = matlab.codetools.requiredFilesAndProducts(fun)';4
list0   = unique([list1;list2;list3;list4]);
listALL = ReadNames(cd);

dnm = true(size(list0));
for i=1:length(dnm)
   [~,~,ext]=fileparts(list0{i});
   if ~strcmp(ext,'.m')
       dnm(i)=false;
   end
end
usedM = list0(dnm);

dnm = true(size(listALL));
for i=1:length(dnm)
   [~,~,ext]=fileparts(listALL{i});
   if ~strcmp(ext,'.m')
       dnm(i)=false;
   end
end
allM = listALL(dnm);

% %%
% 
% for i=1:length(myfiles)
%     cd('D:\Dropbox\SeismicHazard - BETA');
%     mf = myfiles{i};
%     [A,B,C]=fileparts(mf);
%     cd(A)
%     movefile([B,C],'D:\Dropbox\SeismicHazard - BETA\purge')
% end
%%




