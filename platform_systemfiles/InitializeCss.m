function handles=InitializeCss(handles)

handles.Rbin              = createObj('Rbin');
handles.Mbin              = createObj('Mbin');
handles.opt.MagDiscrete   = {'uniform' 0.1};

n = max(handles.sys.branch(:,1:3),[],1);
for i=1:n(3)
    handles.sys.mrr2(i) = process_truncexp  (handles.sys.mrr2(i) , handles.opt.MagDiscrete);
    handles.sys.mrr3(i) = process_truncnorm (handles.sys.mrr3(i) , handles.opt.MagDiscrete);
    handles.sys.mrr4(i) = process_yc1985    (handles.sys.mrr4(i) , handles.opt.MagDiscrete);
end

handles.pop_branch.Value=1;
handles.pop_branch.String=compose('Branch %i',handles.sys.isREG);
handles.pop_site.String=handles.h.id;


fig = handles.figure1;
delete(findall(fig,'type','line'));
delete(findall(fig,'tag','patch'));
delete(findall(fig,'tag','legend_ax1'));
set(findall(fig,'type','axes'),'ColorOrderIndex',1);
set(fig, 'WindowButtonMotionFcn', '');
drawnow
handles.xlabel1  = xlabel(handles.ax1,'Sa(T)','fontsize',8);
handles.ylabel1  = ylabel(handles.ax1,'\lambda Sa(T)','fontsize',8);
handles.param    = {2,'1','0.01 - 2','0.04 - 0.4','5.5 - 8',1,'0 - 1000','350 - 1350','5 - 150','0 - 3000','0 - 1000','32','1000','700 - 0.01','100'};
handles.Tcss     = [0.010;0.050;0.075;0.100;0.150;0.200;0.250;0.300;0.400;0.500;0.750;1.000;1.500;2.000];
handles.AEP      = logsp(1,1e-6,10)';
handles.corrV    = 7;
handles.flatfile = '';
handles.Npreselect = 0;
handles.im1        = [];
handles.lambda1    = [];
handles.lambda2    = [];
handles.T          = [];
handles.eq         = [];
