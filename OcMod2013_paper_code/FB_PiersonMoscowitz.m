function [S,Hs] = FB_PiersonMoscowitz(omega,Tm,moment_no);

if nargin==2
   moment_no = 0;
end

g        = 9.81;
Hs       = g.*(0.4.*Tm./(2.*pi)).^2;
aa       = 8.1e-3*g^2;    % Ochi (1998)
bb       = 5/4;
T        = 2*pi./omega;
S        = aa.*omega.^(moment_no-5).*exp(-bb.*(T./Tm).^4);
