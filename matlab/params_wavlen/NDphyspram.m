function y=NDphyspram(n);
%% THIS PROGRAM STORES ALL THE PHYSICAL PARAMETERS NEEDED FOR THE PROBLEM.
%% CALL: y=NDphyspram(n);
%%
%% IF n==[] => y={h,H_dim,th,[L_ice,T_ice,m_ice]}
%% which are the default values for h, the thickness of the ice that we are
%% nondimensionalising wrt, H_dim is the default water depth,
%% th is the angle of incidence of the incoming wave, L_ice and T_ice are the
%% characteristic length and time of the ice we are nondimensionalising wrt,
%% and since the nondimensional parameter mu=m_ice/lam^(1/4), where 
%% m_ice=rho_ice*h/rho_water/L_ice,
%% so it's sometimes handy to have m_ice around.
%%
%% IF n~=[], & Y=[E g rho_wtr rho_ice nu],
%% THEN n=0 =>y=Y, else y=Y(n).


E        = 5e9;%% Pa
g        = 9.81;%% m/s^2
rho      = 1025;%% kg/m^3
rho_ice  = 922.5;%% kg/m^3
nu       = .3;

if isempty(n)==1
  h      = 1;
  T      = 5;
  H_dim  = 1000;
  th     = 0;%% set default values for these parameters.
  D      = E*h^3/12/(1-nu^2); 
  L_ice  = (D/rho/g)^.25; 
  T_ice  = sqrt(L_ice/g);
  m_ice  = rho_ice*h/rho/L_ice;
  y      = {h,T,H_dim,th,[L_ice,T_ice,m_ice]};
else
  y0  = [E g rho rho_ice nu];
  if n==0
    y = y0;
  else
    y = y0(n);
  end
end
