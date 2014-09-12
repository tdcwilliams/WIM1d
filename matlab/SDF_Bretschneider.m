function S  = SDF_Bretschneider(omega,sdf_prams); 
%% CALL: S  = SDF_Bretschneider(omega,sdf_prams); 
%% sdf_prams={peak period, significant wave height, moment number};

Tm    = sdf_prams{1};
om_m  = 2*pi/Tm;
Hs    = sdf_prams{2};

if length(sdf_prams)==2
   moment_no   = 0;
else
   moment_no   = sdf_prams{3};
end

T  = 2*pi./omega;
f1 = 5/16*Hs^2*om_m^4;
f2 = omega.^(moment_no-5);
f3 = exp(-1.25*(T./Tm).^4);
S  = f1*f2.*f3;
