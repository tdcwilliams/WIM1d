function [Lmiz,Dmiz,output,x_miz]=...
   wim1d_ideal_V20120808(WIM,Tm,h,prof,AC_option,prams_in,xtra_in)
%% CALL: [Dmiz,Lmiz,output]=...
%%    wim1d_ideal_V20120713(WIM,Tm,h,prof,AC_option,...
%%       prams_in,xtra_in);
%% - possibility of using split power law distribution (FSD_CHG==1);
%%    otherwise use Dany''s small floe fractal model;
%% - removed stress condition;
%% - variable advection/attenuation scheme
%%    (see definition of 'adv_method' variable below
%%     adv_method==5 is best - Lax-Wendroff with Superbee flux limiting),
%%    but waves must go through the broken ice that they produce
%%    (for all CFL);
%%   - for checking purposes there is the possibility of not requiring this
%%      (DO_INIT==0)(in input 'prams_in')
%%      (only makes sense with CFL==1 though -
%%       with CFL<1 it is approximately DO_INIT==1 anyway);
%% - fixed attenuation model:
%%    the scattering attenuation is from the A1 model with Y*=5.49GPa 
%%    (Vernon's formula with vbf=[brine vol frac]=0.1)
%%    the routine it uses is:
%%       ALPfxn_E549_cheb2.m;
%%
%%    AC_option = [visc_rp,USE_Wsq],
%%       where visc_rp=[eddy viscosity parameter]=\Gamma;
%%    *damping comes from:
%%       RPget_lam_dmpg.m;
%%    *USE_Wsq==0, W=1; USE_Wsq==1, W=k_ice/k_wtr;
%%       USE_Wsq==2, W=|T|k_ice/k_wtr;
%%     The factor W determines the map S->S*W^2 which gets the ice spectrum
%%       from the water spectrum S;
%% prams_in  = [FSD_CHG,VS_strainc,SAME_SPEED,ICE_LAMB,CFL,DO_INIT]
%% - can use Vernon's breaking strain estimate also (VS_strainc==1);
%%    this also removes fatigue from stress condition;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% extra inputs:
%% xtra_in  = {dx,strain_fac,D_init};
%% strain_c = strain_fac*strainc_default (depends on value of VS_strainc);
%% unless they are specified, dx = 5000m and strain_fac = 1;

dx          = 5000;
strain_fac  = 1;
D_init      = 500;

if nargin==0%%use example inputs;
   disp(['CALL: [Dmiz,Lmiz,output]=',...
            'wim1d_ideal_V20120713(WIM,Tm,h,prof,AC_option,'...
               'prams_in,xtra_in)']);
   WIM         = 3;
   Hs_inc      = 3;
   Tm          = [9.5 Hs_inc];
   c           = .75;
   h           = [3 c];
   prof        = 'exp';
   if 0%%test 'realistic' mode of h,c input;
      X_end    = 450e3;%m
      X_edge   = 45e3;%m
      %%
      [xgrid,hice,cice,edge]  = get_xhc(X_edge,X_end,dx,h(1),c,prof);
      h  = {hice,cice};
   end
   visc_rp     = 10;
   USE_Wsq     = 1
   AC_option   = [visc_rp USE_Wsq];
   %%
%  dx          = 5e3;
%  strain_fac  = 1;
%  D_init      = 500;
   xtra_in     = {dx,strain_fac,D_init};
   %%
   DO_ATTEN_2nd   = 1;
   FSD_CHG        = 0;
   VS_strainc     = 1;
   SAME_SPEED     = 1;
   ICE_LAMB       = 1;
   CFL            = 1;
   DO_INIT        = 1;
   prams_in       = [FSD_CHG,VS_strainc,SAME_SPEED,ICE_LAMB,CFL,DO_INIT]
end

if exist('xtra_in')
   if ~isempty(xtra_in{1})
      dx = xtra_in{1};
   end
   if length(xtra_in)>=2
      if ~isempty(xtra_in{2})
         strain_fac  = xtra_in{2};
      end
   end
   if length(xtra_in)>=3
      if ~isempty(xtra_in{3})
         D_init   = xtra_in{3};
      end
   end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('WIM')
   WIM   = 3;%% ISM method;
             %% WIM==1, WGM, stress & strain criteria;
             %% WIM==2, WGM, strain criterion only;
end

if ~exist('Tm')
   Tm = 7;%% peak wave period
end

if ~exist('h')
   h  = 2;%% ice thickness
end

%% Ice concentration
if iscell(h)
   hcd   = h;
else
   c  = 0.75;
   if length(h)==2
      c  = h(2);
      h  = h(1);
   end
   hcd   = [];
end

if ~exist('prof')
   prof  = [];
end
if isempty(prof)
%  prof        = 'const';%% constant ice thickness
   prof='exp';%% thickness increases away from ice edge
end

visc_rp  = 10;
USE_Wsq  = 0;
if ~exist('AC_option')
   AC_option   = [visc_rp USE_Wsq];
else
   visc_rp  = AC_option(1);%%eddy viscosity
   if length(AC_option)==2
      USE_Wsq  = AC_option(2);
   end
end

%% prams  = [FSD_CHG,VS_strainc,SAME_SPEED,ICE_LAMB,CFL,DO_INIT]
prams = [0 1 1 1 1 1];%%default parameters
if exist('prams_in')
   prams(1:length(prams_in))  = prams_in;
end
FSD_CHG     = prams(1);%% Floe size distribution. 1: do use split power law floe size distribution
                       %% 0: use small floe fractal distribution;
VS_strainc  = prams(2);%% Breaking strain model. 1: do use Vernon's estimate for the breaking strain;
                       %% 0: use Dany's model;
SAME_SPEED  = prams(3);%% 1: waves travel at same speed - no dispersion; 0: waves disperse;
ICE_LAMB    = prams(4);%% 1: correct wavelength in ice; 0: water wavelength;
CFL         = prams(5);%% SAME_SPEED=1: a=CFL*dx/dt; SAME_SPEED=1: dt=CFL*dx/max(a);
                       %% CFL=1 gives no damping by the numerical advection scheme of wim1d_ideal_spec;
%FIX_DT      = prams(6);%% 1: sets dt as a constant and changes speed to fit CFL;
                       %% 0: sets speeds to their usual values and changes dt to fit CFL;
DO_INIT  = prams(6);
%dt_fixed    = 260;
%%
if 1
   USE_DL   = 2;%% Use similar criterion to Vaughan & Squire (2011)
                %%  or Langhorne et al (2001);
   P_crit   = exp(-1);
   %P_crit   = exp(-2);
elseif 1
   USE_DL   = 1;%% use Dave Leslie's (Stats, Uni of Bristol)
                %%  suggestion about the breaking criterion;
   P_crit   = 0.5;
else
   USE_DL   = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Advection scheme:
adv_method  = 5;%%Lax-Wendroff/Superbee
i_start     = 4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if length(Tm)==1%%use Pierson-Moscowitz spectrum;
   [dum,Hs] = FB_PiersonMoscowitz(1,Tm,0);
else%% - else use Bretschneider;
   Hs       = Tm(2);
   Tm(2)    = [];
end

DO_PLOT  = 0;
DO_REP   = (nargin==0);
%edge     = 10;%%ice starts at i=10;
i_tst    = -1;
%i_tst    = 20
%i_tst    = 76

if 0%%check against wim1d_ideal_spec_advexact.m
   SAME_SPEED  = 1;
   i_tst       = 0;
   disp('warning: not using input SAME_SPEED');
end

%% also give progress report every 'reps' time
%%  steps;
reps  = round(10/CFL);
t0    = 24*60^2*rem(now,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%maxNumCompThreads(4);
%format short g 
%format compact

if DO_REP
   disp('Initialization')
end

%% Environment
% homedir = '/home/nersc/dany';
% wavedir = [homedir,'/waves'];
% mizdir  = [homedir,'/validation/miz'];

%% Grid%
conc_crit   = 1e-1;

if isempty(hcd)
   if c<conc_crit
      Dmiz     = NaN;
      Lmiz     = NaN;
      output   = NaN;
      x_miz    = NaN;
      return;
   end
   X_end    = 450e3;%m
   X_edge   = 45e3;%m
   %%
   [xgrid,hice,cice,edge]  = get_xhc(X_edge,X_end,dx,h,c,prof);
   nx                      = length(xgrid);
   %%
   Dmax  = D_init*ceil(cice);
else
   hice  = hcd{1};
   cice  = hcd{2};
   if length(hcd)==3
      Dmax  = hcd{3};
   else
      Dmax  = D_init*ceil(cice);
   end
   nx    = length(hice);
   xgrid = dx*(0:nx-1)';
   %%
   edge     = find(cice>conc_crit);
   if isempty(edge)
      Dmiz     = NaN;
      Lmiz     = NaN;
      output   = {0*xgrid,xgrid,[]};
      x_miz    = NaN;
      return;
   else
      edge     = edge(1);
      X_edge   = xgrid(edge);
   end
end

%% Waves;
fmin  = .042;
fmax  = 1/2.5;
%%
om1   = 2*pi*fmin;
om2   = 2*pi*fmax;
nw    = 31;% NB needs to be odd for Simpson's rule;
dom   = (om2-om1)/(nw-1);
om    = om1+(0:nw-1)'*dom;
T     = 2*pi./om;

%% Get integration weights for each freq;
wt_simp            = 2+0*T;
wt_simp([1 end])   = 1;
wt_simp(2:2:end-1) = 4;
wt_int             = dom/3*wt_simp;

%% Ice thickness
dh      = 0.1;
hh      = 0.4:dh:3.7;
nh      = length(hh); %#ok<NASGU>

%% Model Parameters
rhoice   = 922.5;        % Ice density           [kg/m³]
rho      = 1025.0;       % Water density         [kg/m³]
rhoave   = 0.5.*(rho+rhoice);
g        = 9.81;         % Gravity               [m/s²]

young    = 5e9;              % Young's modulus       [Pa]
poisson  = .3;             % Poisson's ratio
lame_lam = young*poisson/(1+poisson)/(1-2*poisson);
lame_mu  = young/2/(1+poisson);
                        % Lame constants \lambda & \mu  [Pa]
flex_rig_coeff = young/12/(1-poisson^2);
     % flexural rigidity = flex_rig_coeff*(thickness)^3
     % [J/m^3; NB J=kg m^2/s^2]

if VS_strainc
   vbf      = .1;
else
   salt     = 5;            % Ice salinity          [psu]
   temp     = -10;          % Ice temperature       [oC]

   % Brine volume fraction (Frankenstein and Gardner 1967)
   vbf      = 1e-3*salt.*(49.185./abs(temp) + 0.532);
end

% Flexural strength (Timco and O'Brien 1994)
sigma_c = 1.76e6.*exp(-5.88.*sqrt(vbf)); % [Pa]

if VS_strainc
   mu       = 1;
   young    = 1e9*(10*(1-3.51*vbf)-1);       % [Pa]
   poisson  = .3;
else
   % Fatigue, for use in stress condition (Langhorne et al. 1998)
   mu       = 0.6;            %                       [-]
   %%
   young    = 5e9;            % Young's modulus       [Pa]
   poisson  = .3;             % Poisson's ratio
end
Dchg_coeff  = (pi^4/48*young/12/(1-poisson^2)/rho/g)^.25;
   %% D_c=(h^3*pi^4/48*Y/12/rho/g/(1-nu^2))^.25
   %% no flexural breaking below this length

if ~VS_strainc%  Endurance limit (Langhorne et al. 1998)
   strain_c = 3.5e-5;       %                       [-] 
else
   %strain_c = sigma_c/young/(1-poisson^2);
   strain_c = sigma_c/young;
end
strain_c = strain_fac*strain_c;

% Parameters for the floe size distribution
Dmin        = 20;           %                       [m]
xi          = 2;            %                       [-]
fragility   = 0.9;          %                       [-]



% disp('min wave speed = ')
% disp(min(ag))
% disp('max wave speed = ')
% disp(max(ag))

%CFL   = max(ag)*dt/dx;
%if ( CFL>1 )
%    disp('***  Violation CFL   ***')
%    disp(['*** Reduce time step by factor of ',num2str(max(ag)*dt/dx),' ***'])
%    return
%end

%% Ice conditions
% Idealized
%prof     = 'exp';
%h        = 1;


%% USE ICE WAVELENGTH &/OR SPEED;
%% ALSO PRECOMPUTE ATTEN COEFF'S;
ag_ice   = zeros(nw,nx);
wlng_ice = zeros(nw,nx);
alp_all  = zeros(nw,nx);
dmpg_all = zeros(nw,nx);

%% Water wavelength and wave speed as a function only of wave period;
wlng     = g.*T.^2./(2.*pi);
ap       = sqrt(g.*wlng./(2.*pi));     % Phase speed
ag       = ap./2;                      % Group speed
agmax    = max(ag);

if SAME_SPEED
   %dt    = dt_fixed;
   %ag(:) = CFL*dx/dt;%%CFL=1 for testing;
   ag(:) = agmax;
   dt    = CFL*dx/agmax;
   if ICE_LAMB==2
      ICE_LAMB = 1;%%just use ice wavelengths not speeds;
   end
end

%% calc ice wavelengths and speeds (if ICE_LAMB==2);
pramsRP  = [young,poisson,rho,rhoice,g];%%for Robinson-Palmer dispersion relation;

for i=1:nx
   
   [lam_ice,dmpg_all(:,i),ag_ice]   = RPget_lam_dmpg(hice(i),om,pramsRP,visc_rp);

   if ICE_LAMB==2
      wlng_ice(:,i)  = lam_ice;
      %[ag_ice,wlng_ice(:,i)]  = GEN_get_ice_groupvel(hice(i),2*pi./om);
      ag_eff(:,i)             = cice(i)*ag_ice+(1-cice(i))*ag;
      agmax                   = max(agmax,max(ag_eff(:,i)));
   elseif ICE_LAMB==1
      ag_eff(:,i)    = ag;
      wlng_ice(:,i)  = lam_ice;
      %wlng_ice(:,i)  = GEN_get_ice_wavelength(hice(i),2*pi./om);
   else
      ag_eff(:,i)    = ag;
      wlng_ice(:,i)  = wlng;
   end
   %[alp_all(:,i),dmpg_all(:,i)]  = ALPfxn_allACoptions(om,hice(i),AC_option);
   if USE_Wsq==0
      alp_all(:,i)   = ALPfxn_E549_cheb2(om,hice(i));
      Wsq(:,i)       = 1+0*om;
   elseif USE_Wsq==1
      alp_all(:,i)   = ALPfxn_E549_cheb2(om,hice(i));
      Wsq(:,i)       = (wlng./lam_ice).^2;
   elseif USE_Wsq==2
      [alp_all(:,i),Tsq]   = ALPfxn_E549_cheb2(om,hice(i));
      Wsq(:,i)             = Tsq.*(wlng./lam_ice).^2;
   end
end

%[T,wlng_ice(:,edge)/2],pause

if ~SAME_SPEED
%   if FIX_DT%%fix dt and scale speeds;
%      dt       = dt_fixed;
%      ag_eff   = CFL*dx/dt*ag_eff/agmax;
%
%%     ageff_tst   = ag_eff(:,edge-2:edge+1),pause
%%     ag_max_tst  = [T,max(ag_eff,[],2)],pause
%%     tstep_wi    = CFL*dx./ag_eff(:,edge-2:edge+1),pause
%%     tstep_w     = CFL*dx./max(ag_eff,[],2),pause
%   else
      dt = CFL*dx/agmax;
%   end
end
%nsteps  = 250;
nsteps  = 1000;
time    = 0:dt:nsteps*dt;
nt      = length(time);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% freq & spatially dependent (since the speed can depend on space also) time step;
t_glob         = 0;%%global time passed;
dt_freq        = 0*ag_eff;
t_freq         = dt_freq;%%freq dep amount of time passed;
%%
dt_freq(:,1)   = dx./max(ag_eff,[],2);%% freq dep time step
                                   %% (take max speed to get min time step);
for i=2:nx
   dt_freq(:,i)   = dt_freq(:,1);
end
%cfl_tst     = [T,ag_eff(:,edge-1).*dt_freq(:,edge-1)/dx,...
%                ag_eff(:,edge).*dt_freq(:,edge)/dx],pause

%% NB this is indep of space, but because advection is done as a loop over i,
%% we need to have it as a matrix over i
%% NB 2. advection done all at once in hycom so only a vector there
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Memory allocation
%A                 = zeros(nt,nx,nw);
S_atten           = zeros(nx,nw);
Ac                = zeros(nx,nw);

%% Incident wave: Bretschneider spectrum
moment_no   = 0;
%S           = FB_Bretschneider(2*pi./T,Tm,Hs,moment_no); 
S           = SDF_Bretschneider(2*pi./T,{Tm,Hs,moment_no});%, pause;


if 0%% test we have enough integration points:
%   om'
%   wt_int'
   %S'
   variance_wtr   = wt_int'*S;
   %[sum(wt_int), om(end)-om(1)]
   tst_Hs         = [Hs,4*sqrt(variance_wtr)]
   pause;
end

% Boundary conditions (steady state sea: i=1 is same until t=Ntopups*dt)
%if 1%%original setup;
%   Ntopups     = nt;
%   Nwavecells  = 1;
%else
%   Ntopups     = 15;
%   Nwavecells  = edge-1;
%end
%for i=1:Nwavecells
%   S_atten(i,:)   = S;
%end
for i=1:i_start-1
   S_atten(i,:)   = S;
end

%B     = A;
%S_adv = S_atten;

if DO_REP
   % Display parameters
   disp('------------------------------------')
   disp(['WIM        = ',num2str(WIM)]);
   disp(' ');
   %%
   disp(['sigma_c    = ' num2str(sigma_c) ' Pa']);
   disp(['fragility  = ' num2str(fragility)]);
   disp(['fatigue    = ' num2str(mu)]);
   disp(['strain_c   = ' num2str(strain_c)]);
   if iscell(h)
      disp(['average h  = ' num2str(mean(hice(edge:end))) ' m ' prof]);
      disp(['conc       = ' num2str(mean(cice(edge:end)))]);
   else
      disp(['h          = ' num2str(h) ' m ' prof]);
      disp(['conc       = ' num2str(c)]);
   end
   disp(['Tm         = ' num2str(Tm) ' s']);
   disp(['Hs         = ' num2str(Hs) ' m']);
   disp(['AC_option  = ' num2str(AC_option)]);
   disp(' ');
   %%
   disp(['FSD_CHG    = ' num2str(FSD_CHG)]);
   disp(['VS_strainc = ' num2str(VS_strainc)]);
   disp(['SAME_SPEED = ' num2str(SAME_SPEED)]);
   disp(['ICE_LAMB   = ' num2str(ICE_LAMB)]);
   disp(['CFL        = ' num2str(CFL)]);
%   disp(['FIX_DT     = ' num2str(FIX_DT)]);
%  disp(['Ntopups    = ' num2str(Ntopups)]);
%  disp(['Nwavecells = ' num2str(Nwavecells)]);
   disp('');

   %% Integration
   disp('Integrating ...')
end

Lmiz0 = -1;
Dmiz0 = -1;

Hs_i              = zeros(nx,1);
Hs_i(1:i_start-1) = Hs;
n_get_Hs          = 10:5:40;
output3           = cell(1,length(n_get_Hs)+1);
output3{1}        = dt*n_get_Hs;
%%
Dave     = Dmax;
alp_dim  = 0*alp_all;
%            alp_dim(w,i)   = alpha;

%_dummy  = S_atten;
%%NB S_dummy(1,:) never gets overwritten so no need to top up waves;

for n = 1:nt
   %% 2. do advection;
   for w=1:nw
      S_atten(:,w)   = ADV_outer(S_atten(:,w),ag_eff(w,i),dx,dt,adv_method);
   end


%  t_glob                  = t_glob+dt;
%  S_dummy(i_start:end,:)  = S_atten(i_start:end,:);%%SDF from last time step;
%  S_dummy(1:7,:)
%  S_atten(1:7,:)
%  n,pause
   for i = i_start:nx-1
      %% need to compute 4 integrals over period;
      variance_stress = 0;
      variance_strain = 0;
      mom0            = 0;
      mom2            = 0;
      mom4            = 0;
      % ---------------------------------------------------------
      % 0. define no of floes you'd expect to meet as you travel
      % in a line.
      %cice1d   = 0;
      %if cice(i)>0
      %   cice1d         = cice(i)/Dave(i);
      %   if n==1
      %      alp_dim(:,i)   = alp_all(:,i)*cice1d+2*cice(i)*dmpg_all(:,i);
      %      %% from now on, this is only updated after breaking;
      %   end
      %end%% check if ice is present;


      BREAK_CRIT  = 0;
      for w=1:nw%% NB important to loop in order of increasing frequency
                %% (decreasing wavelength) with WGM method;


         % ----------------------------------------------------------
         % 1. Update wavelength and wave speed
         % if it depends on ice conditions.
         % Other wise precompute before the time loop.
      
         k_ice    = 2*pi/wlng_ice(w,i);
         W        = sqrt(Wsq(w,i));

         % ------------------------------------------------------
         % 2. Compute critical wave height for floe breaking
         %    (only used for WGM method)
         if (n == 1)&(WIM<3)
            Ac1   = wlng_ice(w,i).^2.*strain_c./...
                     (2*W.*pi.^2.*hice(i));
                     %% strain criterion
            if 1
               Ac2_fac  = 2/(1-poisson^2);
            else
               Ac2_fac  = 1;
            end
            Ac2   = Ac2_fac*2.*pi.*hice(i).^2.*mu.*...
                     sigma_c./(3*W.*g.*rhoave.*wlng_ice(w,i).^2);
                     %% stress criterion
                     %% - NB factor of 2/(1-nu^2) difference from Dumont et al (2011);
                     %% - plates are stiffer than beams and consider only a half-wavelength
                     %%   for the 3pt loaded beam analogy;
            %%
            if WIM==1%%both strain & stress;
               Ac(i,w)  = min(Ac1,Ac2);
            else%%only stress
               Ac(i,w)  = Ac1;
            end
         end

         % ----------------------------------------------------------
         % 3. Wave advection
%         A(n,i,w) = A(n-1,i,w) + ...
%             ag_eff*dt*(A(n-1,i-1,w) - A(n-1,i,w))/dx;
%         B(n,i,w) = B(n-1,i,w) + ...
%             ag_eff*dt*(B(n-1,i-1,w) - B(n-1,i,w))/dx;
         %%
%        ds       = 0;%%determines the amount of attenuation;
%        ds_glob  = ag_eff(w,i)*(t_glob-t_freq(w,i));
         %%
%        if i==4&w==10
%           n
%           tt = [t_freq(w,i),t_glob]
%           pause
%        end

         %% ADVECT EACH FREQ UNTIL t_freq(w,i)>=t_glob;
         %% THEN INTERP VALUE BACK TO t_glob (GLOBAL/BREAKING TIME);
%         if t_freq(w,i)<t_glob%%??should this be <=??
%                               %% NO - we want to stop if
%                               %% t_freq = t_glob (perfect) or if
%                               %% t_freq > t_glob (need to then interp back);
%            S_dummy(i,w)   = S_atten(i,w);
%            ds             = ds+ag_eff(w,i)*dt_freq(w,i);
%            S_atten(i,w)   = S_dummy(i,w) + ...
%                              ds/dx*(S_dummy(i-1,w) - S_dummy(i,w));
%%           if i>=2 & i<=4 & w==10
%%              old   = S_dummy(i-1:i,w)
%%              new   = S_atten(i-1:i,w)
%%              [n i t_freq(w)],pause;
%%           end
%            t_freq(w,i)      = t_freq(w,i)+dt_freq(w,i);
%
%            % ----------------------------------------------------------
%            % 2. Update energy attenuation coefficient
%            % from a look-up table
%            % Convert from energy to AMPLITUDE by dividing by 2;
%            % convert to dimensional coefficient;
%            alpha          = alp_all(w,i)/2*cice1d+dmpg_all(w,i)*cice(i);
%            alp_dim(w,i)   = alpha;
%         else
%            % don't update attenuate coefficient
%            % as it will stop waves getting out of cell for low CFL;
%            alpha    = alp_dim(w,i);
%            ds_glob  = ds_glob+dx;
%         end
%
%%         S_adv(n,i,w) = S_adv(n-1,i,w) + ...
%%             ag_eff*dt*(S_adv(n-1,i-1,w) - S_adv(n-1,i,w))/dx;
%         t0       = t_freq(w,i)-dt_freq(w,i);
%         S0       = S_dummy(i,w);
%         t1       = t_freq(w,i);
%         S1       = S_atten(i,w);
%         S_glob   = S0 + (S1-S0)/(t1-t0)*(t_glob-t0);
%
%         % ----------------------------------------------------------
%         % 3. Wave advection & attenuation
%%        A(n,i,w)       = A(n,i,w)*exp(-alpha*ag_eff*dt);
%         S_atten(i,w)   = S_atten(i,w)*exp(-2*alpha*ds);
%         S_glob         = S_glob*exp(-2*alpha*ds_glob);
         if DO_INIT==1
            if cice(i-1)==0
               alpha = 0;
            else
               cice1d   = cice(i-1)/Dave(i-1);
               alpha    = alp_all(w,i-1)*cice1d+2*dmpg_all(w,i-1)*cice(i-1);
            end
            Atten = alpha*ag_eff(w,i-1)*dt;
            %% NB Breaking has already been done in cell i-1 so 'Atten'
            %% is the attenuation due to the broken ice (FOR ALL CFL);
            CC    = ag_eff(w,i-1)*dt/dx;
         else
            if cice(i)==0
               alpha = 0;
            else
               cice1d   = cice(i)/Dave(i);
               alpha    = alp_all(w,i)*cice1d+2*dmpg_all(w,i)*cice(i);
            end
            Atten = alpha*ag_eff(w,i)*dt;
            %% NB Breaking has not yet been done in cell i so 'Atten'
            %% is the attenuation due to the ice ahead 
            %% (BUT ONLY FOR CFL **VERY** CLOSE TO 1);
            CC    = ag_eff(w,i)*dt/dx;
         end
         S_atten(i,w)   = S_atten(i,w)*exp(-Atten);
         %if 0
         %   S0             = S_dummy(i-1,w)*exp(-Atten);
         %   S1             = S_dummy(i,w);
         %   S_atten(i,w)   = CC*S0+(1-CC)*(S1-S0);
         %else
         %   dS             = feval(adv_fxn,S_dummy(i+adv_jvec,w),CC);
         %   S_atten(i,w)   = (S_dummy(i,w)+dS)*exp(-Atten);
         %end
%        if (i==n+1)&(n>8)&(w==10)
%           disp('tstS')
%           [n i]
%           {alp_all(w,i-1),Dave(i-1),...
%              alp_all(w,i-1)*cice(i-1)/Dave(i-1),2*dmpg_all(w,i-1)*cice(i-1)}
%           Atten
%           S0,S1,[S_dummy(i-1,w),S_atten(i,w)],pause
%        end
         S_glob         = S_atten(i,w);
         
%        if ((n==3)|(n==4))&(i==4)*(w==10)
%           n
%           ttt   = [t0,t1,t_glob]
%           bet   = (t_glob-t0)/(t1-t0)
%           SSS = [S_dummy(i,w),S_atten(i,w),S_glob]
%           pause
%        end


%        if i==10&w==10
%           n,w,i,om(w),S_atten(i,w),W
%           [wlng_ice(i,w),wlng(w)],pause
%        end
         if 0%i==11 & n>=15
            wt_int
            tstSint  = [S_atten(i,w),i,n,T(w)]
            pause
         end
         %% energy spectrum of wtr particle at i-th node;

         % ----------------------------------------------------------
         % 5. Compute strain, & integrate;
         %    also need m_2
         %    (Cartwright & Longuet-Higgins, 1956);
         mom0 = mom0 + wt_int(w)*S_glob*W^2;
         if hice(i)>0
%             stress_density    = S_atten(i,w)*W^2*...
%                                  (3*g*rhoave*wlng_ice(w)^2)^2/...
%                                   (2*pi*hice(i)^2)^2;
             strain_density    = S_glob*W^2*...
                                  (k_ice^2*hice(i)/2)^2;
             SD_tst(w,:)   = [strain_density W k_ice^2/2 hice(i)];
%             variance_stress   = variance_stress +...
%                                    wt_int(w)*stress_density;
             variance_strain   = variance_strain +...
                                    wt_int(w)*strain_density;
             %%
             mom2    = mom2 +...
                       + wt_int(w)*S_glob*W^2*om(w)^2; 

            %%
            if WIM<3
               %% DO BREAKING INSIDE FREQ (w) LOOP;
               wamp           = sqrt(2*om(w)*S_glob);
               BREAK_CRITwgm  = ( (wamp>Ac(i,w))&(wlng_ice(w,i)/2)>Dmin );
                  %% is amplitude large enough to break floes
                  %% & will new Dmax be >Dmin?
               BREAK_CRITwgm  = BREAK_CRITwgm&( wlng_ice(w,i)/2<Dmax(i) );
                  %% Dmax needs to decrease for breaking to occur;
               if FSD_CHG
                  %% breaking can only occur if Dmax isn't too small already
                  %% (then it won't break as it's a rigid body effectively);
                  Dchg           = Dchg_coeff*hice(i)^.75;
                  BREAK_CRITwgm  = BREAK_CRITwgm&( Dchg<Dmax(i) );
               end
               if BREAK_CRITwgm
                  Dmax(i)     = wlng_ice(w,i)/2;
                  BREAK_CRIT  = 1;
               end%%
               if i==i_tst&(w==3|w==4)&mom0>0
                  wa_tst   = [i w wamp Ac(i,w) BREAK_CRITwgm],pause
               end
            end%% end WIM<3 test
         end%%end test for ice present
         % ----------------------------------------------------------
      end%%end w loop (omega)
      Hs_i(i) = 4*sqrt(mom0);
      do_tst  = (i==i_tst);
      %do_tst  = (i==(n+1))&(n>=8);

      if hice(i)>0 & mom2>0
         Nwaves   = dt/2/pi*sqrt(mom2/mom0);
         T_crit   = 2*pi*sqrt(mom0/mom2);
         %%
         sig_strain  = 2*sqrt(variance_strain);
         %%
%         Pstress  =...
%            FB_wave_height_Raleigh_truncdist(sigma_c,sig_stress);
         Pstrain  =...
            FB_wave_height_Raleigh_truncdist(strain_c,2*sig_strain);

         if WIM==3
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% FLOE BREAKING FOR WIM==3:
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if FSD_CHG
               Dchg        = Dchg_coeff*hice(i)^.75;
               BREAK_CRIT  = (Dchg<Dmax(i));%% is current Dmax<D_c? (can't break anymore)
            else
               BREAK_CRIT  = 1;
            end

            if USE_DL==1
               %% use Dave Leslie's (Stats, Uni of Bristol) suggestion
               %% as other criterion may be too tough;
               P_nobreak   = exp(Nwaves*log(1-Pstrain));%%(1-Pstrain)^Nwaves
               BREAK_CRIT  = BREAK_CRIT&( (1-P_nobreak)>P_crit );
            elseif USE_DL==2
               %% Use similar criterion to Vaughan & Squire (2011)
               %%  or Langhorne et al (2001);
               BREAK_CRIT  = BREAK_CRIT&( Pstrain>P_crit );
            else
               P_crit      = 1/Nwaves;
               BREAK_CRIT  = BREAK_CRIT&( Pstrain>=P_crit );
            end

            if BREAK_CRIT
               %% use crest period to work out wavelength
               %% - biggest floes will break into half this;
               if ICE_LAMB==0
                  wlng_crest = g.*T_crit.^2./(2.*pi);
               else
                  wlng_crest =...
                     RPget_lam_dmpg(hice(i),2*pi/T_crit,pramsRP,0);
                     %old wavelength used incorrect
                     %Young's modulus:
                     %GEN_get_ice_wavelength(hice(i),T_crit);

               end

               Dc = wlng_crest/2;
               if Dc >= Dmin & Dmax(i)>Dc
                  Dmax(i) = Dc;
               else
                  BREAK_CRIT  = 0;%% don't need to change Dmean
                                  %% as we haven't changed Dmax;
               end
            end
         end%%end WIM==3 breaking;

         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %% NOW WORK OUT AVERAGE FLOE SIZE IF BREAKING HAS OCCURRED:
         if BREAK_CRIT
            if FSD_CHG
               gam1        = 1.15;
               gam2        = 2.5;
               Dchg        = Dchg_coeff*hice(i)^.75;
               FSD_prams   = {gam1,gam2,Dchg};
            else
               FSD_prams   = {fragility,xi};
            end%%end FSD choice
            Dave(i)  = FSD_Dave(Dmin,Dmax(i),FSD_prams,FSD_CHG);
            %dv = [i,Dave(i)],pause
         end%%end breaking
%'hey',do_tst,pause 
         %%
         if do_tst
%            subplot(1,2,1);
%            plot(T,squeeze(A(n,i,:))), hold on;
%            plot(T,squeeze(B(n,i,:)),'k'), hold off;
%            set(gca,'yscale','log');
%            %%
%            subplot(1,2,2);
%            plot(T,log10(squeeze(S_atten(i,:)))), hold on;
%            plot(T,log10(squeeze(S_adv(n,i,:))),'k'),hold off;
%            ylim([-10 0]);
            disp(['Hs=',num2str(4*sqrt(mom0)),'m, Tw=',num2str(T_crit),'s']);
            %disp([num2str(4*sqrt(mom0)),'m, ',num2str(T_crit),'s']);%,mom4
            %crest_period,Dc
            %stresses = [sigma_c,sig_stress,sigma_c<sig_stress]
            if WIM==3
               Nwaves
               strains  = [strain_c,sig_strain,sig_strain/strain_c]
               variance_strain
               probs    = [Pstrain,1-(1-Pstrain)^Nwaves,P_crit]
               probs2   = [Pstrain,exp(-2*strain_c^2/sig_strain^2)]
               %SD_tst   = S_atten([i-1 i],:)'
               %Hs_tst   = [4*sqrt(SD_tst'*wt_int),4*sqrt(SD_tst'*(wt_int.*(wlng./wlng_ice(:,i)).^2))]%,SD_tst
            end
            disp('n, t, i, BC, Dmax')
            if BREAK_CRIT, disp('BREAKING'), end;
            disp([n,n*dt/60,i,BREAK_CRIT,Dmax(i)]);%,Dc
            disp('pause'),pause
         end%%'do_tst' check
      end%% check for ice and waves present;
   end%%end spatial loop i;
   
   %% calc Dmiz & Lmiz;
   dd    = Dmax(:)';
   jmiz  = find(dd>0&dd<D_init);
   Lmiz  = dx/1e3*length(jmiz);
   Dmiz  = max(dd(jmiz));

   if i_tst>0&Lmiz>0
      jmiz2    = edge-1:(edge-1+round(1e3*Lmiz/dx)+1);
      jmiz20   = 1:(edge-1+round(1e3*Lmiz/dx)+1);
      %Hs0   = Hs_i(1:i_tst+4)
      Hs0   = [jmiz20',Hs_i(jmiz20)]

      %Dtst  = [Dave(jmiz2),Dmax(n-1,jmiz2)',Dmax(n,jmiz2)']
      %Dtst     = [jmiz2',Dave(jmiz2),Dmax(jmiz2)]
      nDLmiz   = [n,Dmiz,Lmiz]
      pause
   elseif i_tst>0
      Hs0      = [(1:i_tst+4)',Hs_i(1:i_tst+4)]
      nLwaves  = [n, length(find(Hs0))+1]
      pause
   end

   if mod(n,reps)==0%% progress report;
      %%
      t1    = 24*60^2*rem(now,1);
      if DO_REP
         disp([n nt]);
         %disp([num2str(t1-t0),' mins taken']);
         %%
         Dmax_edge   = dd(10:20)'
         disp(['Dmiz, Lmiz:']);
         disp([Dmiz,Lmiz]);
         %{num2str(Tm),num2str(Hs),num2str(AC_option)}
         disp(['Tm = ',num2str(Tm),...
               ', Hs = ',num2str(Hs),...
               ', AC_option = ',num2str(AC_option)]);
         disp(['Time taken (s): ',num2str(t1-t0)])
      end
      %%
      if (ICE_LAMB<2)&(abs(1-Lmiz/Lmiz0)<1e-6&abs(1-Dmiz/Dmiz0)<1e-6)
         break;
      end
      %%
      Lmiz0 = Lmiz;
      Dmiz0 = Dmiz;
   end%% end progress report;


   n_ = find(n==n_get_Hs);
   if ~isempty(n_)
      output3{1+n_}   = Hs_i;
   end

end%% end time loop;

if ~exist('Dmiz','var')
   Dmiz  = NaN;
elseif isempty(Dmiz)
   Dmiz  = NaN;
end

%% Termination
% ----------------------------------------------------------------------

if DO_REP
   % ----------------------------------------------------------------------
   % Display results
   disp('Displaying results :')
   disp(['Peak period      = ',num2str(Tm),' s']);
   disp(['Sig wave height  = ',num2str(Hs),' m']);
   disp(['Mean thickness   = ',...
            num2str(mean(hice(edge:end))),' m']);
   disp(['Lmiz             = ' num2str(Lmiz) ' km'])
   disp(['Dmiz             = ' num2str(Dmiz) ' m'])
   disp(['AC_option        = ' num2str(AC_option)])

   disp(['SAME_SPEED       = ' num2str(SAME_SPEED)])
   disp(['profile shape    = ' prof])
   disp(['USE_DL           = ' num2str(USE_DL)])
   disp(['log[P_crit]      = ' num2str(log(P_crit))])

   disp('------------------------------------')
end

output   = {Dmax,xgrid,output3};
x_miz    = X_edge/1e3+Lmiz;

if isempty('xedge')
   x_miz = NaN;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xgrid,hice,cice,edge]  = get_xhc(X_edge,X_end,dx,h,c,prof);

xgrid    = (0:dx:X_end)';
nx       = length(xgrid);
edge     = 1+ceil(X_edge/dx);
%%
hice           = xgrid;
cice           = xgrid;
hice(1:edge-1) = 0;
cice(1:edge-1) = 0;
cice(edge:end) = c;
%%
if strcmp(prof,'const')
   hice(edge:end) = h;
elseif strcmp(prof,'exp')
   L_h   = 60e3;%m
   h0    = .1*h;
   dh    = h-h0;
   %%
   ex_h           = exp(-(xgrid(edge:end) - xgrid(edge))./L_h);
   hice(edge:end) = h0+dh.*(1 - ex_h);
   if 0
      plot(xgrid,hice);
      pause;
   end
else
   prof
   disp('Wrong profile: exp or const')
   return
end
