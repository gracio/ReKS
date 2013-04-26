function dLogHalf = perturbBeta(beta0,s,U,sqrtD)
%
% put a ceiling on half-life.
% watch out for S = 1
% sqrtD : a column vector

  % eq 7.21
  tmp = log(2)/(s * log(s) * (beta0 * log(s) - log(2))); 
  
  q = sqrtD.^(-1) .* U; 
  
  % eq 7.28
  dLogHalf = tmp * (2 * q * q' - ... 
		   s*q.^2 * ones(1, length(q)) - ... 
		   s*(ones(length(q),1) *(q.^2)')); 
 