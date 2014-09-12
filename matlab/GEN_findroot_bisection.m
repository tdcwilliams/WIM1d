function [x,exitflag]=...
GEN_findroot_bisection(fxn,interval,varargin)

tol=1e-12;
%%
%z=varargin{:}%%extra arguments for function to be zeroed.
i0=interval(1);
i1=interval(2);

fl=feval(fxn,i0,varargin{:});
fh=feval(fxn,i1,varargin{:});

%% check there is a root in interval
if (fl > 0 & fh > 0) | (fl < 0 & fh < 0 )
  x=[];
  exitflag=0;
  disp('warning (GEN_findroot_bisection.m):');
  disp('no root in interval');
  return;
end

%% if either endpoint is "close enough", just use that for root
if (fl == 0)
  x = i0;
  exitflag=1;
  return;
end
if (fh == 0)
  x = i1;
  exitflag=2;
  return;
end

%% define search direction to be in dirn of increasing p
if (fl < 0.0)
  xl=i0;
  xh=i1;
else
  xl=i1;
  xh=i0;
end

x=0.5*(xl+xh);
dx=(xh-xl)/2;
%  {x,dx,xl,xh}

Nits=0;
critter=1;

while critter
  Nits=Nits+1;
%    if Nits==36;{x,dx,xl,xh},end;
  p=feval(fxn,x,varargin{:});
  if p==0
    exitflag=3;
    return;
  elseif p>0
    xh=x;
  else
    xl=x;
  end
  dx=(xh-xl)/2;
  ERR=abs(dx/x);
  critter=( abs(dx)>tol & ERR>tol );
  x=xl+dx;%{x,dx,xl,xh}
end
exitflag=4;