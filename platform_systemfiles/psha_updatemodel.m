function [sys,opt,h]=psha_updatemodel(matfile,fname)

sys = [];
opt = [];
h   = [];
if nargin==1 && exist([matfile,'.mat'],'file')==2
    load(matfile,'sys','opt','h')
elseif nargin==1 && exist([matfile,'.mat'],'file')~=2
    [sys,opt,h]= loadPSHA([matfile,'.txt']);
%     D=what('platform_defaultmodels'); save([D.path,'\',matfile],'h','sys','opt');
elseif nargin==2
    [sys,opt,h]= loadPSHA(fname);
end