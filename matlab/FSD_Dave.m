function Dave  = FSD_Dave(Dmin,Dmax,FSD_prams,FSD_CHG)
%% CALL: Dave  = FSD_Dave(Dmin,Dmax,FSD_prams,FSD_CHG)
%% OUTPUT:
%% *Dave = average floe size;
%% INPUTS:
%% *Dmin = min floe size possible;
%% *Dmax = max floe size possible (determined by wave field);
%% *if FSD_CHG==0 -> Dumont et al (2011) model;
%%   *FSD_PRAMS = {fragility xi};
%% *if FSD_CHG==1 -> Split power law - like in Toyota et al (2010);
%%   *FSD_prams = {gam1 gam2 Dchg};
%%     *gam1 is exponent for small floes;
%%     *gam2 is exponent for small floes;
%%     *Dchg is D where regime shift (change in exponent) occurs;

D_top = 200;
if Dmax>D_top
   Dave  = Dmax;
   return;
end
%%
if Dmax<Dmin
   Dave  = Dmax;
   %Dave  = Dmin;
   return;
end

if FSD_CHG
%  gam1        = 1.15;
%  gam2        = 2.5;
%  Dchg        = Dchg_coeff*hice(i)^.75;
   gam1        = FSD_prams{1};
   gam2        = FSD_prams{2};
   Dchg        = FSD_prams{3};
   if 0%%truncate the long floe PDF at Dmax;
      P0          = [];%% make PDF continuous;
      PDF_prams   = {Dmin,gam1+1,Dchg,gam2+1,Dmax};
      Dave        = PDF_SplitPowerLaw_trunc_pdf([],PDF_prams); 
   else%%
      PDF_prams   = {Dmin,gam1+1,Dchg,gam2+1};
      Pmax        = 0.95;
      %Dmax(n-1,i)
      PDF_prams   = PDF_SplitPowerLaw_getP0(Dmax,Pmax,PDF_prams);
      Dave        = PDF_SplitPowerLaw_pdf([],PDF_prams);
   end
else%% Dumont et al (2011) model;
%  fragility   = 0.9;
%  xi          = 2;
   fragility   = FSD_prams{1};
   xi          = FSD_prams{2};
   Dave        = subfn_floe_scaling(fragility,xi,Dmin,Dmax);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Dave = subfn_floe_scaling(f,xi,Dmin,Dmax)
% This function computes the average floe size within a grid cell as a
% function of the maximum floe size using a bounded fractal renormalization
% group method.

% We suggest to use Dmin >= 20. Below that value, there is no scattering
% by the floes for periods larger than 6 s (it's probably viscous however).

%%% LGB CHANGE / Oct 2013

if 1
 gam  = 2 + log(f)/log(xi);
 Dave = gam*(Dmax*((Dmin/Dmax)^gam)-Dmin)/(1-gam);
else
 M  = floor(log2(Dmax/Dmin));
 if isfinite(M) && M > 0
    m = 0:M;
    N = (1 - f).*(xi.^2.*f).^m;
    ND = N.*Dmax./(xi.^m);
    Dave = sum(ND)./sum(N);
 else
    Dave = Dmin;
 end
end

Dave = max(Dave,Dmin);
