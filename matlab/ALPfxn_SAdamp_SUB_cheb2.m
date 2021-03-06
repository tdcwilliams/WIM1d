function [alp,damping]=ALPfxn_SAdamp_SUB_cheb2(om,h);
%% made from 'attenLuke4Mar2011LongFA_cheb2.dat' 
%% & 'dampingLuke4Mar2011LongFA_cheb2.dat'
%% which were obtained by running ALPfxn_LB4Mar2011_cheb.m

if nargin==0
   om = 2*pi/10;
   h  = 1;
end


if h==0
   alp      = 0*om;
   damping  = alp;
   return;
end


Ncheb_f  = 10;%% deg of poly in period/freq; 
Ncheb_h  = 10;%% deg of poly in thickness;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INPUT AND OUTPUT FILES:

%% copy cheb coeff's from
%% 'attenTWrslPL_c90_cheb2.dat'
%% into function 'alp_coeffs_fn'
%% at the bottom of this file;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SET INPUT TYPES:

%% X_INPUT  = 1 for period;
%% X_INPUT  = 2 for freq;
%% X_INPUT  = 3 for omega;
X_INPUT  = 2;

%% AC_INPUT  = 1 if input AC is energy AC;
%% AC_INPUT  = 2 if input AC is amplitude AC;
AC_INPUT = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RANGE OF VALIDITY VARIES FROM RUN TO RUN,
%%  SO NEED TO CHECK THIS EACH TIME;
hmax  = 5;
hmin  = .1;
xmin  = .042;  %% Tmax=25s
xmax  = .4;    %% Tmin=2.5s
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EVERYTHING FROM HERE ON SHOULD BE THE SAME
%%  FOR EACH RUN;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

zfac  = 1;
jout  = find(h>hmax | h<hmin);
if ~isempty(jout)
   if 0
      disp(sprintf('h outside range [%d,%d]',hmin,hmax));
      return;
   elseif h>hmax
      h  = hmax;
   else
      zfac  = h/hmin;
      h     = hmin;
   end
end
%%

if X_INPUT==1
   x     = (2*pi)./om;
   xstr  = 'period';
elseif X_INPUT==2
   x     = om/(2*pi);
   xstr  = 'freq';
else
   x  = om;
   xstr  = 'omega';
end

jout  = find(x>xmax | x<xmin);
if ~isempty(jout)
   x(jout)
   disp([xstr,...
     sprintf(' outside range [%d,%d]',xmin,xmax)]);
   return;
end

alp_coeffs  = alp_coeffs_fn(Ncheb_f,Ncheb_h); 
damp_coeffs = damp_coeffs_fn(Ncheb_f,Ncheb_h); 
Limits      = [xmin,xmax,hmin,hmax];

%% DO INTERPOLATION:
y_coeffs    = alp_coeffs+1i*damp_coeffs;
y           = OP_chebinterp2d(x,h,y_coeffs,Limits);

%% want to OUTPUT ENERGY attenuation coefficients:
alp         = zfac*AC_INPUT*exp(real(y));
damping     = zfac*AC_INPUT*exp(imag(y));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = alp_coeffs_fn(nfc,nhc);

out     = zeros(nfc+1,nhc+1);

out(:)  =...
[ -1.8501579e+00	
   4.9101877e+00	
  -1.9944048e+00	
   7.3571385e-01	
  -2.2388319e-01	
   4.3700266e-02	
   1.4351128e-02	
  -2.5879087e-02	
   1.9640453e-02	
  -1.0363481e-02	
   3.2663770e-03	
   2.7695669e+00	
  -1.2063711e+00	
  -5.0784397e-02	
   3.7480698e-01	
  -2.6425166e-01	
   9.1308132e-02	
   6.9547408e-03	
  -4.3420410e-02	
   4.6092399e-02	
  -3.1540644e-02	
   1.5837220e-02	
  -1.2036662e+00	
   3.3680013e-01	
   1.9808329e-01	
  -2.3814800e-01	
   8.0319397e-02	
   3.7840528e-02	
  -6.0775037e-02	
   3.5031731e-02	
  -6.0416843e-03	
  -8.9097036e-03	
   1.5726785e-02	
   6.4898367e-01	
  -1.2109648e-01	
  -1.5806140e-01	
   1.3688927e-01	
  -3.3480367e-03	
  -6.1277527e-02	
   4.1973575e-02	
  -5.6549417e-03	
  -1.1632590e-02	
   1.2089542e-02	
  -6.4837898e-03	
  -3.7936273e-01	
   5.3998050e-02	
   1.1371674e-01	
  -7.5373235e-02	
  -2.6600942e-02	
   5.5102220e-02	
  -2.1563646e-02	
  -8.0610915e-03	
   1.2680156e-02	
  -5.9534129e-03	
  -1.1215486e-03	
   2.3485146e-01	
  -2.8291759e-02	
  -8.0908767e-02	
   3.8793499e-02	
   3.4616020e-02	
  -3.9923118e-02	
   4.7538859e-03	
   1.4474075e-02	
  -1.0138135e-02	
   1.0790529e-03	
   3.2584851e-03	
  -1.5053499e-01	
   1.2525344e-02	
   5.7500633e-02	
  -1.8200037e-02	
  -3.1466412e-02	
   2.5868808e-02	
   4.0116291e-03	
  -1.4768318e-02	
   6.5896646e-03	
   1.8758681e-03	
  -3.7027599e-03	
   1.0190258e-01	
  -2.2429425e-04	
  -4.0608910e-02	
   6.9273313e-03	
   2.5656968e-02	
  -1.5656594e-02	
  -7.7986690e-03	
   1.2440697e-02	
  -3.0118289e-03	
  -3.7434000e-03	
   3.4392031e-03	
  -6.7257320e-02	
  -4.5018646e-03	
   2.6230343e-02	
  -1.4532063e-03	
  -1.8832335e-02	
   8.5065359e-03	
   8.4602431e-03	
  -9.2579927e-03	
   4.6304848e-04	
   4.2940931e-03	
  -2.7030121e-03	
   3.8848120e-02	
   3.8380170e-03	
  -1.4649663e-02	
  -2.4250313e-04	
   1.1576818e-02	
  -3.9017295e-03	
  -6.5477149e-03	
   5.7063587e-03	
   7.0878167e-04	
  -3.4892451e-03	
   1.7178745e-03	
  -2.3909477e-02	
  -2.8499561e-03	
   8.8697690e-03	
   5.1143991e-04	
  -7.5240865e-03	
   1.6398133e-03	
   5.0596886e-03	
  -3.4272682e-03	
  -1.3717067e-03	
   2.8594988e-03	
  -9.5172114e-04  ];

%% factor of 2 was missing in original definition of alp
%% (ie ALPfxn_LB_LFA4Mar2011_cheb.m used to be amp AC not energy AC)
%% so added log(2) to constant term;
%% NB this IS NOW done in ALPfxn_LB_LFA4Mar2011_cheb.m;
%out(1)   = out(1)+log(2);


function out = damp_coeffs_fn(nfc,nhc);

out      = zeros(nfc+1,nhc+1);
out(:)   = ...
[ -1.2233019e+01	
   1.4519969e-01	
  -3.0148111e-01	
   2.3742724e-01	
  -1.5423088e-01	
   9.6480359e-02	
  -5.9199747e-02	
   3.4991675e-02	
  -1.9572012e-02	
   9.8009109e-03	
  -4.2839433e-03	
  -1.5073001e+00	
  -7.3553740e-01	
   5.2598861e-01	
  -2.5163648e-01	
   7.1908038e-02	
   6.2572098e-03	
  -2.8200328e-02	
   2.8164817e-02	
  -2.0418112e-02	
   1.1440691e-02	
  -6.1785019e-03	
   4.9484433e-01	
   4.0946407e-01	
  -2.2304337e-01	
   3.0487214e-02	
   5.7320175e-02	
  -5.9445032e-02	
   3.4480403e-02	
  -1.1908350e-02	
  -1.3558438e-03	
   5.6946671e-03	
  -6.1196024e-03	
  -2.3318404e-01	
  -2.3101814e-01	
   8.7637783e-02	
   3.3088361e-02	
  -5.9667711e-02	
   3.1089041e-02	
  -1.3657105e-03	
  -1.1877202e-02	
   1.2965055e-02	
  -8.8384312e-03	
   4.1215163e-03	
   1.2414779e-01	
   1.3891949e-01	
  -3.5294341e-02	
  -3.8149307e-02	
   3.6841709e-02	
  -5.7353116e-03	
  -1.3700209e-02	
   1.4190742e-02	
  -6.9953807e-03	
   1.2273073e-03	
   2.4705757e-03	
  -7.1257123e-02	
  -8.8753470e-02	
   1.3123439e-02	
   3.0938423e-02	
  -1.8761548e-02	
  -6.3752909e-03	
   1.4583399e-02	
  -7.6341361e-03	
  -7.1560475e-04	
   3.7549782e-03	
  -3.8705399e-03	
   4.2180028e-02	
   5.7105361e-02	
  -3.5010876e-03	
  -2.2144227e-02	
   8.3109209e-03	
   9.4134848e-03	
  -1.0648099e-02	
   1.9406586e-03	
   4.2310236e-03	
  -4.4944972e-03	
   2.2127542e-03	
  -2.5778203e-02	
  -3.7181215e-02	
  -7.3173365e-04	
   1.5259532e-02	
  -2.8733395e-03	
  -8.8122776e-03	
   6.4546850e-03	
   1.3593037e-03	
  -4.7253396e-03	
   3.1414945e-03	
  -8.4790829e-05	
   1.5835451e-02	
   2.3867144e-02	
   2.2131525e-03	
  -1.0122262e-02	
   3.6525172e-04	
   6.8649923e-03	
  -3.3519888e-03	
  -2.5467776e-03	
   3.7813389e-03	
  -1.5179202e-03	
  -1.1687535e-03	
  -8.8847762e-03	
  -1.3757294e-02	
  -2.0636176e-03	
   5.9459058e-03	
   4.5093754e-04	
  -4.3833673e-03	
   1.4296488e-03	
   2.2285197e-03	
  -2.3517064e-03	
   4.3926964e-04	
   1.3171873e-03	
   5.3308274e-03	
   8.4306204e-03	
   1.7888482e-03	
  -3.6963192e-03	
  -7.1646906e-04	
   2.8757952e-03	
  -4.4809105e-04	
  -1.8383124e-03	
   1.4040431e-03	
   1.7212564e-04	
  -1.2364047e-03  ];
