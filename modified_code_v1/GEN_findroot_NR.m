function [x,unfinished_roots]=...
  GEN_findroot_NR(fxn,guess,varargin)
%% CALL:[x,exitflag]=GEN_findroot_NR(fxn,guess,varargin)
%% x is the root of 'fxn' nearest to 'guess'
%% NB 'guess' can be a vector
%%
%% varargin contains any extra arguments needed
%%  to calculate 'fxn' (eg. fxn=@sin)
%%
%% unfinished_roots=[] => all have converged
%% unfinished_roots~=[] => some haven't converged -
%%  the ones in the unfinished_roots vector.

tol               = 1e-9;
MAXITS            = 300;%% need this to stop infinite loops
unfinished_roots  = (1:length(guess))';
%%
j  = 0;
x0 = guess;

while ~isempty(unfinished_roots) & j<MAXITS
  [p,dp] = feval(fxn,x0,varargin{:});
  dx     = -p./dp;
  x      = x0+dx;%pause
  %%
  j   = j+1;
  x0  = x;
  %%
  critter            = (abs(dx)>tol)|(abs(p)>1);
  unfinished_roots   = find(critter);
end

if j==MAXITS & ~isempty(unfinished_roots)
  disp('warning (GEN_findroot_NR.m): root not converged');
  x(unfinished_roots)   = NaN;
end
