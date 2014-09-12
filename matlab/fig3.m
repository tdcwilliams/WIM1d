%clear;
%%1x2 plot: eg of advection
if 1
   CFL_vec     = [[1 .99 .95 .9  .7] ];
   %CFL_vec     = [[1 .99 .95 .9  .7 .5] ];
end

if 1
   visc_rp  = 13;
   dir0     = 'out2';
else
   visc_rp  = 10;
   dir0     = 'out';
end
if USE_Wsq==2 
   figname  = [dir0 '/fig3B_Tsq.eps'];
   filename = [dir0 '/fig3B_Tsq.mat'];
else
   figname  = [dir0 '/fig3B.eps'];
   filename = [dir0 '/fig3B.mat'];
end
opt      = 2;
%%
Hs    = 4;
Tm    = [9.5 Hs];
c     = .8;
h     = [2 c];
prof  = 'exp';
%%
dx       = 5e3;
D_init   = 500;

if exist(filename)
   load(filename);
else

   if 1%~exist('h_VS')%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      disp('VS runs');
      dx_VS    = 5e3;%.1e3;
      DO_2nd   = 0;
      DO_INIT  = [1 0];
      %%
      for j=1:length(DO_INIT)
         xtra_in              = {dx_VS,1,D_init};
         prams_in             = [DO_2nd DO_INIT(j) 0];
         [Lmiz,Dmiz,hcd,Hs]   = Lmiz_analytic(Tm,h,prof,xtra_in,prams_in);
         %Lmiz,Dmiz
         %%
         Lmiz_VS(j)  = Lmiz;
         Dmiz_VS(j)  = Dmiz;
         %%
         Dmax_VS{j}     = hcd{3};
         xgrid_VS       = hcd{4}/1e3;%% [km]]
         h_VS           = hcd{1};
         Hs_VS{j}       = Hs;
         staggerVS(j)   = (0==DO_INIT(j))*dx_VS/1e3;
      end
   end%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   Lmiz_VS,Dmiz_VS
   
   if opt==1%[[N1][N2]]
      DO_INIT_vec = [[1 ], [0  0  0   0  0 ] ];
      CFL_vec     = [[1 ], [1 .99 .95 .9 .8] ];
   elseif opt==2%%get rid of N2 results
      DO_INIT_vec = [[0  0  0   0   0  0] ];
      CFL_vec     = [[1 .99 .95 .9  .7 .5] ];
   else
      DO_INIT_vec = [[1 1],  [0  0   0  0 0]  ];
      CFL_vec     = [[1 .1], [1 .99 .95 .9 .8] ];
   end
   if 1%~exist('Dmax_AAS')%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      disp('AAS runs');
      staggerAAS  = (0==DO_INIT_vec)*dx/1e3;
      naas        = length(DO_INIT_vec);
      %%
      WIM         = 3;
      %AC_option   = 2;
      USE_Wsq     = 2;
      AC_option   = [visc_rp USE_Wsq];
      prof        = 'exp';
      %%
      FSD_CHG     = 0;
      VS_strainc  = 1;
      SAME_SPEED  = 1;
      ICE_LAMB    = 1;
      for j=1:naas
         xtra_in   = {dx,1,D_init};
         prams_in  = [FSD_CHG,VS_strainc,SAME_SPEED,ICE_LAMB,CFL_vec(j),DO_INIT_vec(j)];
         %[Lmiz,Dmiz,output]   = wim1d_ideal_V20120621(WIM,Tm,h,prof,AC_option,prams_in,xtra_in);
         %[Lmiz,Dmiz,output]   = wim1d_ideal_V20120713(WIM,Tm,h,prof,AC_option,prams_in,xtra_in);
         [Lmiz,Dmiz,output]   = wim1d_ideal_V20120808(WIM,Tm,h,prof,AC_option,prams_in,xtra_in);
         %%
         %Lmiz,Dmiz
         Lmiz_AAS(j) = Lmiz;
         Dmiz_AAS(j) = Dmiz;
         %%
         Dmax_AAS{j} = output{1};
         xgrid_AAS   = output{2}/1e3;
         Hs_AAS{j}   = output{3};
         %Hs_VS{j}    = Hs;
      end
   end%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   save(filename,'Lmiz_AAS','Dmiz_AAS','Dmax_AAS','staggerAAS','xgrid_AAS','Hs_AAS',...
           'Lmiz_VS','Dmiz_VS','Dmax_VS','staggerVS','xgrid_VS','Hs_VS','h_VS',...
              'CFL_vec');
end
%Lmiz_AAS,Dmiz_AAS
if opt==1
   colin = {'^k','.k','--b','-b','ob','*b'};
elseif opt==2
   %colin = {'.k','--b','-b','ob','*b','^k','vg'};
   ls    = {'o','-','-','-','-','-','-'};
   lc    = {[0 0 0],.4*[0 1 1],.75*[0 1 1],[0 0 .5],[0 0 1]};
   lc    = lc([1 2 4 3 5]);
else
   colin = {'^k','vk','.k','--k','-k','ob','.b'};
end

str   = {'Hs, m','Dmax, m'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax    = .17;
X     = [0 250];
ay    = .08;
Y     = [0 500];
Y1    = [0 2];
ay34  = .12;
Y34   = [5e-2 5];
%%
lw = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(1,2,2);
%%
x1 = xgrid_VS+staggerVS(1);
x2 = xgrid_VS+staggerVS(2);
y1 = Hs_VS{1};
y2 = Hs_VS{2};

PQ(1)    = plot(x1,y1,'color',[1 0 0],'linewidth',lw,'linestyle','-');
leg2{1}  = 'A0';
hold on;
PQ(2)    = plot(x2,y2,'color',[.5 0 0],'linewidth',lw,'linestyle','-');
leg2{2}  = 'A1';
%plot(x1,y1,'r',x2,y2,'--r');
%%
set(gca,'yscale','log');
hold on;
%%
np = 5;
%t1 = output{3}{1}(np)/60^2
if opt==0
   j1 = 1;
   j2 = 3;
elseif opt==1
   j1 = 1;
   j2 = 2;
else
   j1       = 4;
   CFL_figb = CFL_vec(j1)
   j2       = 1;
end
t1 = Hs_AAS{j1}{1}(np)/60^2
x1 = xgrid_AAS+staggerAAS(j1);
x2 = xgrid_AAS+staggerAAS(j2);
y1 = Hs_AAS{j1}{np+1};
y2 = Hs_AAS{j2}{np+1};
%%
ssc   = 'iC=';
PQ(3)    = plot(x2,y2,'color',[0 0 0],'linewidth',1,'linestyle','o');
leg2{3}  = [ssc '1'];
PQ(4)    = plot(x1,y1,'color',[0 0 0],'linewidth',4.5,'linestyle','.');
%leg2{4}  = ['C = ' num2str(CFL_figb)];
leg2{4}  = [ssc num2str(CFL_figb)];
%plot(x1,y1,'.k',x2,y2,'ok');
legend(PQ,leg2,'location','northeast');
GEN_proc_fig('x, km',str{1});
%%
x     = X(1)+ax*(X(2)-X(1));
Y0    = log10(Y34);
y     = 10^( Y0(1)+ay34*(Y0(2)-Y0(1)) );
txt   = text(x,y,'(b)');
GEN_font(txt);
hold off;
axis([X Y34]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%subplot(2,2,4);
%x1 = xgrid_VS+staggerVS(1);
%x2 = xgrid_VS+staggerVS(2);
%plot(x1,Hs_VS{1},'r',x2,Hs_VS{2},'--r');
%set(gca,'yscale','log');
%hold on;
%%
%np = 6;
%%t2 = output{3}{1}(np)/60^2
%t2 = Hs_AAS{j1}{1}(np)/60^2
%x1 = xgrid_AAS+staggerAAS(j1);
%x2 = xgrid_AAS+staggerAAS(j2);
%y1 = Hs_AAS{j1}{np+1};
%y2 = Hs_AAS{j2}{np+1};
%%%
%plot(x1,y1,'.k',x2,y2,'ok');
%GEN_proc_fig('x, km',str{1});
%%%
%axis([X Y34]);
%x     = X(1)+ax*(X(2)-X(1));
%Y0    = log10(Y34);
%y     = 10^( Y0(1)+ay34*(Y0(2)-Y0(1)) );
%txt   = text(x,y,'(d)');
%GEN_font(txt);
%hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax    = .73;
X     = [0 250];
ay    = .12;
Y     = [0 500];
Y1    = [0 2];
ay34  = .2;
Y34   = [5e-2 5];
%%
lw = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(1,2,1);
PP(1)    = plot(xgrid_VS,Dmax_VS{1},'color',[1 0 0],'linewidth',lw,'linestyle','-');
leg{1}   = 'A0';
hold on;
PP(2)    = plot(xgrid_VS,Dmax_VS{2},'color',[.5 0 0],'linewidth',lw,'linestyle','-');
leg{2}   = 'A1';

for j=2:length(CFL_vec)
   PP(j+2)  = plot(xgrid_AAS,Dmax_AAS{j},'color',lc{j},'linestyle',ls{j},'linewidth',lw);
   %PP(j+2)  = plot(xgrid_AAS,Dmax_AAS{j},colin{j},'linewidth',lw);
   leg{j+2} = [ssc,num2str(CFL_vec(j))];
end
for j=1
   PP(j+2)  = plot(xgrid_AAS,Dmax_AAS{j},'color',lc{j},'linestyle',ls{j},'linewidth',1);
   %PP(j+2)  = plot(xgrid_AAS,Dmax_AAS{j},colin{j},'linewidth',lw);
   leg{j+2} = [ssc,num2str(CFL_vec(j))];
end
%legend(PP,leg,'best')
legend(PP,leg,'location','northwest')
GEN_proc_fig('x, km',str{2});
%%
axis([X Y]);
x     = X(1)+ax*(X(2)-X(1));
Y0    = Y;
y     = Y0(1)+ay*(Y0(2)-Y0(1));
txt   = text(x,y,'(a)');
GEN_font(txt);
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%subplot(2,2,1);
%plot(xgrid_VS,h_VS,'k');
%hold on;
%GEN_proc_fig('x, km','h, m');
%%%
%%axis([X Y1]);
%xlim(X)
%x     = X(1)+ax*(X(2)-X(1));
%Y0    = Y1;
%y     = Y0(1)+ay*(Y0(2)-Y0(1));
%txt   = text(x,y,'(a)');
%GEN_font(txt);
%hold off;

saveas(gcf,figname,'epsc');
%!gv out/fig3.eps &
%if DO_INIT==0
%   !gv out/lendistQA_cnvgnc_DI0.eps &
%else
%   !gv out/lendistQA_cnvgnc_DI1.eps &
%end
