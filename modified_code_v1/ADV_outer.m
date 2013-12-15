function u = ADV_outer(u,c,h,dt,method)
%% CALL: u = ADV_outer(u,c,h,dt,method)
%% 1d advection code;
%% INPUTS:
%% u is a vector of the thing to be advected;
%% c is a scalar speed;
%% h is the spatial resolution;
%% dt is the temporal resolution;
%% method is an integer that specifies the advection scheme:
%%  method == 1%% upwind, 1st order;
%%  method == 2%%Lax-Friedrichs - central diference, 1st order;
%%  method == 3%%Lax-Wendroff - direct space-time, 2nd order;
%%  method == 4%%Lax-Wendroff with von Leer flux limiting;
%%  method == 5%%Lax-Wendroff with Superbee flux limiting;
%% OUTPUT:
%% u is the vector of the thing to be advected - after advection;

if c<0
   %% some of the methods crash for c<0,
   %% so reverse grid, run algorithm, and reverse again;
   u0 = flipud(u);
   u0 = ADV_outer(u0,-c,h,dt,method);
   u  = flipud(u0);
   return;
end


if nargin<5
   method   = 5;
end
f  = flux(u,c,h,dt,method);
u  = u-dt*diffl(f)/h;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = flux(u,c,h,dt,method)

switch method
case 1%% upwind, 1st order
   f  = c*u;
case 2%%Lax-Friedrichs - central diference, 1st order;
   f  = c*sumr(u)/2-h/2/dt*diffr(u);
case 3%%Lax-Wendroff - direct space-time, 2nd order;
   f  = c*sumr(u)/2-c^2*dt/2/h*diffr(u);
case 4%%Lax-Wendroff with von Leer flux limiting;
   theta = diffl(u)./(diffr(u)+3e-14);
   phi   = limiter(theta,1);
   f     = c*u+c/2*(1-c*dt/h)*diffr(u).*phi;
case 5%%Lax-Wendroff with Superbee flux limiting;
   theta = diffl(u)./(diffr(u)+3e-14);
   phi   = limiter(theta,2);
   f     = c*u+c/2*(1-c*dt/h)*diffr(u).*phi;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = diffl(x)
y = [0;diff(x)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = diffr(x)
y = [diff(x);0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = sumr(x)
y = [x(1:end-1)+x(2:end);2*x(end)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function phi = limiter(r,method)
% flux limiter
if nargin<1, method = 2; end
switch method
case 1
   phi = (abs(r)+r)./(1+abs(r));                        % van Leer
case 2
   phi = max(0,max(min(1,2*r),min(r,2)));               % Superbee
end
