function GEN_proc_fig(x_lbl,y_lbl,fontsize,figsize)
%% CALL: GEN_proc_fig(x_lbl,y_lbl,fontsize)
%% x_lbl and y_lbl are strings containing text for the x&y labels;
%% fontsize is the size of the numbers and axes labels;
%% figsize={marg,wid,hei}, which if inputed
%% sets eps (or any other format) fig produced to have margins
%% marg, width wid, and height hei (all in cm);

if ~exist('fontsize')
   fontsize = 17;
end

hold on;
set(gca,'FontName','Times','FontSize',fontsize);
if ~isempty(x_lbl)
   hx = xlabel(x_lbl);
   set(hx,'FontName','Times','FontSize',fontsize,...
         'FontWeight','light');
         %'Interpreter','latex');
end
if ~isempty(y_lbl)
   hy = ylabel(y_lbl);
   set(hy,'FontName','Times','FontSize',fontsize,...
         'FontWeight','light');
         %'Interpreter','latex');
end
hold off;

if exist('figsize')
   marg  = figsize{1};
   wid   = figsize{2};
   hei   = figsize{3};
   GEN_setsize_eps(marg,wid,hei);
end
