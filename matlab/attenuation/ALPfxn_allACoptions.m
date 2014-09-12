function [alp,damping] = ALPfxn_allACoptions(om,h,ACoption)

damping  = 0*om;
if h==0
   alp   = 0*om;
   return;
end

if ACoption==1|ACoption==2
   ALPfxn   = @ALPfxn_SAdamp_SUB_cheb2;%% faster version :)
   %% NB need extra argument 'Dave'
elseif ACoption==3
   ALPfxn   = @ALPfxn_SA_NOSUB_cheb2;
elseif ACoption==4
   ALPfxn   = @ALPfxn_TWrslPL_cheb2;%% faster version :)
elseif ACoption==5;
   if  0
      ALPfxn   = @ALPfxn_KM08_cheb2;%% faster version :)
   else
      alp   = alpha_km08(2*pi./om,h);
      return;
   end
end

% ----------------------------------------------------------
% Retrieve energy attenuation coefficient
if ACoption==2
   [alp,damping]  =...
      feval(ALPfxn,om,h);
else
   alp   = feval(ALPfxn,om,h);
end
