function GEN_setsize_eps(marg,wid,hei)
%% CALL: GEN_setsize_eps(marg,wid,hei)
%% sets eps (or any other format) fig produced to have margins
%% marg, width wid, and height hei (all in cm);

marg0 = .15;%%1.5mm - too low if need labels
if ~exist('marg')
   marg  = marg0;
end
if isempty(marg)
   marg  = marg0;
end

wid0  = 16;
if ~exist('wid')
   wid  = wid0;
end
if isempty(wid)
   wid  = wid0;
end

hei0  = 14;
if ~exist('hei')
   hei  = hei0;
end
if isempty(hei)
   hei  = hei0;
end

ppos  = [marg,marg,[wid,hei]];
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',ppos);
