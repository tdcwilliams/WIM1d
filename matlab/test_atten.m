%%test_atten.m

h_vec       = (.1:.1:5)';%%thicknesses to test
T_vec       = (2:1:25)';%%periods to test
AC_option   = 1;

nh = length(h_vec);
np = length(T_vec);

for j=1:nh
   for r=1:np
      h  = h_vec(j);
      T  = T_vec(r);
      om = 2*pi/T;
      alp(j,r) = ALPfxn_allACoptions(om,h,ACoption);
   end
end
