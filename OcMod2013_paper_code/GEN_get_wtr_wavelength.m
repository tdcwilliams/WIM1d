function wavlen=GEN_get_wtr_wavelength(T,H)
%% CALL: wavlen=GEN_get_wtr_wavelength(T,varargin)
%% calc's open water wavelength ('wavlen') corresponding to a given
%% wave period 'T', and water depth (optional argument - if no
%% value is specified infinite water depth is assumed)

om = 2*pi./T;
g  = NDphyspram(2);%% m/s^2

%% infinite depth wavenumber:
k  = om.^2/g;

if nargin==1
   wavlen   = 2*pi./k;
   return;
elseif H==Inf
   wavlen   = 2*pi./k;
   return;
end

%% if want wavelength for finite depth
%% use the infinite depth 'k' as a seed.
alpha  = 1./k;
for j=1:length(k)
  k(j)   = gen_root_wtr(alpha(j),H,k(j));
end

wavlen   = 2*pi./k;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function k  = gen_root_wtr_v2(alpha,H,guess)

lam      = alpha;
R        = max(1/lam,7/H);
interval = [0 1.2*R];
k        = GEN_findroot_bisection(...
            @disprel,interval,lam,H);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function k  = gen_root_wtr(alpha,H,guess)
%% finds the root of the ice dispersion relation nearest to 'guess'.

tol      = 1e-8;
k0       = guess;
dk       = NR_corr_term(k0,alpha,H);
k        = k0-dk;
critter  = 1;

while critter
  k0        = k;
  dk        = NR_corr_term(k0,alpha,H);
  k         = k0-dk;
  critter   = ( abs(dk)>tol & abs(dk/k)>tol );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dk=NR_corr_term(k,alpha,H)
%% dk=f/f_k, where f has the same zeros as of the dispersion function, 
%% is the correction term in the Newton-Rhapson method for finding zeros in f.

x  = 7.5;
if real(k*H)<=x
  f   = alpha*k*sinh(k*H)-cosh(k*H);
  df  = alpha*k*H*cosh(k*H)+(alpha-H)*sinh(k*H);
else
  f   = alpha*k*tanh(k*H)-1;
  df  = alpha*k*H+(alpha-H)*tanh(k*H);
end
dk = f/df;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f=disprel(k,lam,H)
f  = lam*k*tanh(k*H)-1;
