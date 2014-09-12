function [wavlen,char]  = GEN_get_ice_wavelength(h,T,H)
%% CALL: wavlen=GEN_get_ice_wavelength(h,T,varargin)
%% calc's wavelength ('wavlen') corresponding to a given
%% ice thickness 'h', wave period 'T', and water depth
%% (optional argument - if no
%% value is specified infinite water depth is assumed)

if h==0%% if ice is water, use GEN_get_wtr_period.m
  if nargin==2
    wavlen  = GEN_get_wtr_wavelength(T);
  else
    wavlen  = GEN_get_wtr_wavelength(T,H);
  end
  char   = [];
  return;
end

if nargin==2
   H  = Inf;
end

pram     = NDphyspram(0);
E        = pram(1);%% Pa
g        = pram(2);%% m/s^2
rho      = pram(3);%% kg/m^3
rho_ice  = pram(4);%% kg/m^3
nu       = pram(5);

om       = 2*pi./T;
D        = E*h.^3/12/(1-nu^2);
L        = ((D/rho)./om.^2).^.2;
lam      = g./L./om.^2;
mu       = (rho_ice*h/rho)./L;
del      = lam-mu;
wavlen   = 0*del;
%%
L_ice =(D/rho/g).^.25;
char  = {L,L_ice};


%% infinite depth results:
for j=1:length(del)
  r   = roots([1 0 0 0 del(j) -1]);
  k   = r(find(r>0 & imag(r)==0));
  if H~=Inf
    %% if want wavelength for finite depth,
    %% use the infinite depth 'k' as a seed.
    k = gen_root_ice(del(j),H/L(j),k);
  end
  wavlen(j) = 2*pi*L(j)/k;
end

function k=gen_root_ice(del,H,guess)
%% finds the root of the ice dispersion relation nearest to 'guess'.

tol   = 1e-8;
k0    = guess;
dk    = NR_corr_term(k0,del,H);
k     = k0-dk;
while abs(dk) > tol
  k0  = k;
  dk  = NR_corr_term(k0,del,H);
  k   = k0-dk;
end

function dk = NR_corr_term(k,del,H)
%% dk=f/f_k, where f has the same zeros as of the dispersion function, 
%% is the correction term in the Newton-Rhapson method for finding zeros in f.

Lam   = k^4+del;
Lampr = 5*k^4+del;
x     = 7.5;
if real(k*H)<=x
  f   = Lam*k*sinh(k*H)-cosh(k*H);
  df  = Lam*k*H*cosh(k*H)+(Lampr-H)*sinh(k*H);
else
  f   = Lam*k*tanh(k*H)-1;
  df  = Lam*k*H+(Lampr-H)*tanh(k*H);
end
dk = f/df;
