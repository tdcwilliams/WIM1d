function [P,jbad]=FB_wave_height_Raleigh_truncdist(X,H_s)
%% CALL: P=FB_wave_height_Raleigh_truncdist(X,H_s),
%% where H_s is the significant wave height,
%% and P=probability of getting a wave AMPLITUDE bigger than X;

%% No limit (LGB change)

sig_sq   = (H_s/2)^2;
str_sq   = X.^2;
P        = exp(-2*str_sq/sig_sq);

%% truncated Rayleigh dist, so no rogue waves (P(2*X>2*Hs)=0); 
% 
% P  = 0*X;
% jj = find(X<=H_s);
% %%
% sig_sq   = (H_s/4)^2;
% A        = 1 - exp(-H_s.^2/2/sig_sq);
% P(jj)    = ( exp(-X(jj).^2/2/sig_sq) - exp(-H_s.^2/2/sig_sq) )/A;

return