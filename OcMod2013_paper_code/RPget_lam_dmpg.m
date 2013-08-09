function [lam_ice,damping,ag_ice] = RPget_lam_dmpg(h,om_vec,prams,visc_rp)
%% CALL: [lam_ice,damping,ag_ice] = RPget_lam_dmpg(h,om_vec,prams,visc_rp)
%%  prams = [Y,nu,rho_wtr,rho_ice,g];
%%  dispersion relation (INFINITE DEPTH) is:
%%  f=(D/rho/g)*k^5+(1-alp*d-1i*om*visc_rp/rho/g)*k-alp=0,
%%  where D is flexural rigidity, alp=om^2/g, and visc_rp
%%  is the Robinson & Palmer eddy viscosity parameter;
%% Outputs the ice wavelength and group velocity (from k when visc_rp=0)
%%  and the damping from R&P (=imag(k));
%% TODO: finite depth?

do_test  = 0;

if nargin==0
   h        = 0.5;
   fmin     = .042;
   fmax     = 1/2.5;
   om_vec   = 2*pi*linspace(fmin,fmax,31)';
   %%
   Y        = 5.49e9;
   nu       = 0.3;
   rho      = 1025;
   rho_ice  = .9*rho;
   g        = 9.81;
   %%
   prams    = [Y,nu,rho,rho_ice,g];
   visc_rp  = 10;
   do_test  = 1;
end

Y     = prams(1);%%Young's modulus [Pa]
nu    = prams(2);%%Poisson's ratio
rho   = prams(3);%%density of water [kg/m^3]
rho_i = prams(4);%%density of ice [kg/m^3]
g     = prams(5);%%gravity[m/s^2]
D     = Y*h^3/12/(1-nu^2);

nw       = length(om_vec);
gs1      = om_vec(1)^2/g;%%guess for without viscosity;
gs2      = om_vec(1)^2/g;%%guess for with viscosity;
lam_ice  = 0*om_vec;
ag_ice   = 0*om_vec;
damping  = 0*om_vec;

for w=1:nw
   om    = om_vec(w);
   k_ice = GEN_findroot_NR(@disp_fn_rp,gs1,h,om,prams,0);
   k_rp  = GEN_findroot_NR(@disp_fn_rp,gs1,h,om,prams,visc_rp);
   %%
   gs1   = k_ice;
   gs2   = k_rp;
   %%
   lam_ice(w)  = 2*pi/k_ice;
   damping(w)  = imag(k_rp);

   %% GET GROUP VELOCITY:
   ag_ice(w)   = (5*D*k_ice^4+rho*g-rho_i*h*om^2)/2/om/(rho+rho_i*k_ice*h);

%  %% GET INTRINSIC ADMITTANCE:
%  s(w)  = 
end




if do_test==1
   plot(om_vec, damping,'k');
   if 1
      [alp,dmpg]  = ALPfxn_allACoptions(om_vec,h,2);
      hold on;
      plot(om_vec,dmpg,'r');
      om_vec,[damping,dmpg]
      hold off;
   end
   set(gca,'yscale','log');
   GEN_proc_fig('\omega, s^{-1}','\epsilon, m^{-1}');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [f,df] = disp_fn_rp(k,h,om,prams,visc_rp)
%% CALL: [f,df] = disp_fn_rp(prams,h,om,visc_rp)
%% dispersion relation (INFINITE DEPTH) is:
%% f=(D/rho/g)*k^5+(1-alp*d-1i*om*visc_rp/rho/g)*k-alp=0,
%% where D is flexural rigidity, alp=om^2/g, and visc_rp
%% is the Robinson & Palmer eddy viscosity parameter;
%%
%% prams = [Y,nu,rho_wtr,rho_ice,g];
Y     = prams(1);%%Young's modulus [Pa]
nu    = prams(2);%%Poisson's ratio
rho   = prams(3);%%density of water [kg/m^3]
rho_i = prams(4);%%density of ice [kg/m^3]
g     = prams(5);%%gravity[m/s^2]
d     = rho_i/rho*h;%%draft [m]

alp   = om^2/g;
D     = Y*h^3/12/(1-nu^2);
C5    = D/rho/g;
gam   = om*visc_rp/rho/g;
C1    = 1-alp*d-1i*gam;

f  = C5*k.^5+C1*k-alp;
df = 5*C5*k.^4+C1;
