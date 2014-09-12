function [Dave,Dm,Pm] = FSD_Dave(Dmax,Dmin,lam,f0,model)

if ~exist('model')
   %%0: original fragility model
   %%1: modified to have linear drop from f=f_big to 0
   %%    for D\in [D_min,D']
   %%2: modified to give f=(1-cos(D/D'))^2.5*f_big
   model = 0;
end

if ~exist('Dmax')
   Dmax  = 200;%m
end
if ~exist('Dmin')
   Dmin  = 30;%m
end
if ~exist('lam')
   lam   = 300;%m
end
if ~exist('f0')
   f0 = .9;
end
%%
xi = 2;
if 1
   M  = log(Dmax/Dmin)/log(xi);%=log_xi(Dmax/Dmin)
   %[M,log2(Dmax/Dmin)]
   M  = floor(M);
else
   M  = 40;
end

mm = (0:M)';
Dm = Dmax*(1/xi).^mm;

if model==0
   Nm = (1-f0)*(xi^2*f0).^mm;
   %%
   Ntot  = sum(Nm);
   Pm    = Nm/Ntot;
   Dave  = sum(Pm.*Dm);
   return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fm    = frag(Dm,f0,Dmin,lam,model);
Nm    = 1+0*Dm;
Nm(1) = 1;

for m=0:M-1
   Nm(m+2)  = xi^2*fm(m+1)*Nm(m+1);
   Nm(m+1)  = (1-fm(m+1))*Nm(m+1);
end

Ntot  = sum(Nm);
Pm    = Nm/Ntot;
Dave  = sum(Pm.*Dm);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f=frag(D,f0,Dmin,lam,model)

%%large floes break with prob f0;
f  = f0+0*D;

%%small floes can't break;
jj    = find(D<Dmin);
f(jj) = 0;

%%medium floes have length-dependent frag;
jj    = find((D>=Dmin) & (D<=lam));

switch model
case 1
   f(jj) = (D(jj)-Dmin)/(lam-Dmin)*f0;
case 2
   cc    = cos(pi/2*D(jj)/lam);
   f(jj) = (1-cc).^2.5*f0;
end
