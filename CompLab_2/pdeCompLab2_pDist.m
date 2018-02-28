% This script is written and read by pdetool and should NOT be edited.
% There are two recommended alternatives:
 % 1) Export the required variables from pdetool and create a MATLAB script
 %    to perform operations on these.
 % 2) Define the problem completely using a MATLAB script. See
 %    http://www.mathworks.com/help/pde/examples/index.html for examples
 %    of this approach.
function pdemodel
[pde_fig,ax]=pdeinit;
pdetool('appl_cb',1);
set(ax,'DataAspectRatio',[10 10 1]);
set(ax,'PlotBoxAspectRatio',[1 0.60860860860860855 0.60860860860860855]);
set(ax,'XLimMode','auto');
set(ax,'YLim',[-10 10]);
set(ax,'XTickMode','auto');
set(ax,'YTickMode','auto');

% Geometry description:
pderect([-10 10 3 -3],'Rectangle');
pderect([-3.5 -1.5 1 -1],'Square');
pdeellip(2.25,1.25,0.75,0.75,...
0,'Circ2');
pdeellip(2.25,-1.25,0.75,0.75,...
0,'Circ1');
set(findobj(get(pde_fig,'Children'),'Tag','PDEEval'),'String','Rectangle-Square-Circ2-Circ1')

% Boundary conditions:
pdetool('changemode',0)
pdesetbd(16,...
'neu',...
1,...
'0',...
'0')
pdesetbd(15,...
'neu',...
1,...
'0',...
'0')
pdesetbd(14,...
'neu',...
1,...
'0',...
'0')
pdesetbd(13,...
'neu',...
1,...
'0',...
'0')
pdesetbd(12,...
'neu',...
1,...
'0',...
'0')
pdesetbd(11,...
'neu',...
1,...
'0',...
'0')
pdesetbd(10,...
'neu',...
1,...
'0',...
'0')
pdesetbd(9,...
'neu',...
1,...
'0',...
'0')
pdesetbd(8,...
'neu',...
1,...
'0',...
'0')
pdesetbd(7,...
'neu',...
1,...
'0',...
'0')
pdesetbd(6,...
'neu',...
1,...
'0',...
'0')
pdesetbd(5,...
'neu',...
1,...
'0',...
'-25')
pdesetbd(4,...
'neu',...
1,...
'0',...
'0')
pdesetbd(3,...
'neu',...
1,...
'0',...
'0')
pdesetbd(2,...
'neu',...
1,...
'0',...
'0')
pdesetbd(1,...
'dir',...
1,...
'1',...
'0')

% Mesh generation:
setappdata(pde_fig,'Hgrad',1.3);
setappdata(pde_fig,'refinemethod','regular');
setappdata(pde_fig,'jiggle',char('on','mean',''));
setappdata(pde_fig,'MesherVersion','preR2013a');
pdetool('initmesh')
pdetool('refine')
pdetool('refine')
pdetool('refine')
pdetool('refine')

% PDE coefficients:
pdeseteq(1,...
'1.0',...
'0.0',...
'0',...
'1.0',...
'0:10',...
'0.0',...
'0.0',...
'[0 100]')
setappdata(pde_fig,'currparam',...
['1.0';...
'0.0';...
'0  ';...
'1.0'])

% Solve parameters:
setappdata(pde_fig,'solveparam',...
char('0','91392','10','pdeadworst',...
'0.5','longest','0','1E-4','','fixed','Inf'))

% Plotflags and user data strings:
setappdata(pde_fig,'plotflags',[4 1 4 1 3 1 2 1 0 0 0 1 1 1 0 1 0 1]);
setappdata(pde_fig,'colstring','1e6 + 1/2  * (25.^2 - ux.^2 - uy.^2) ');
setappdata(pde_fig,'arrowstring','[ux; uy]');
setappdata(pde_fig,'deformstring','');
setappdata(pde_fig,'heightstring','1e6 + 1/2  * (25.^2 - ux.^2 - uy.^2) ');

% Solve PDE:
pdetool('solve')
