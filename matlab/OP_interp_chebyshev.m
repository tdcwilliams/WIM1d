function [f,hn]=OP_interp_chebyshev(tt,An)

if ~iscell(An)
  NgC    = length(An)-1;
  hn     = pi/2+zeros(NgC+1,1);
  hn(1)  = pi;
  %%
  f   = 0*tt;

  C0  = 1+f;
  f   = f+An(1)*C0;

  if NgC==0
    return;
  else
    C1   = tt;
    f    = f+An(2)*C1;
  end

  for its=2:NgC
    Cn   = 2*tt.*C1-C0;
    f    = f+An(its+1)*Cn;
    C0   = C1;
    C1   = Cn;
  end

else
  NgC    = An{1};
  hn     = pi/2+zeros(NgC+1,1);
  hn(1)  = pi;
  %%
  f   = zeros(length(tt),NgC+1);
  if NgC<0
    return
  end

  C0     = 1+0*tt;
  f(:,1) = C0;

  if NgC==0
    return;
  else
    C1      = tt;
    f(:,2)  = C1;
  end

  for its=2:NgC
    Cn         = 2*tt.*C1-C0;
    f(:,its+1) = Cn;
    C0         = C1;
    C1         = Cn;
  end
end
