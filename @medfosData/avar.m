 function [a, mu] = avar(obj, Z, p, alpha)
  % Z    : Discrete values of RV
  % p    : probabilities
  % a    : Average value at risk (level alpha)
  % mu   : A subgradient of AVAR_alpha at Z
  n = length(Z);
  [mu, a, exitflag] = linprog(-Z',[],[], ...
     ones(1, n), 1, zeros(n,1), p/alpha);
  assert(exitflag == 1, 'numerical problems');
  a = -a;
  return 