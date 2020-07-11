function handles=mGMPEfromfigures(handles,filename)

handles.rad1.Value=1;
handles.rad2.Value=0;

delete(findall(handles.ax1,'type','line'))
ch=get(handles.panel2,'children');
set(ch(handles.text),'Visible','off')
set(ch(handles.edit),'Visible','off','Style','edit');
str = regexp(filename,'\_','split');
str = str{1};
handles.ax1.ColorOrderIndex = 1;
handles.HoldPlot.Value=1;
handles=mGMPEdefault(handles,ch(handles.text),ch(handles.edit));
handles.epsilon = 0;
fun=str2func(['GMMValidation_',str]);
fun(handles,filename);