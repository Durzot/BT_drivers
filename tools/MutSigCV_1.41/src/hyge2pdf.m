function p=hyge2pdf(k,n,k1,n1)
  p = exp(gammaln(n1+2) - gammaln(k1+1) - gammaln(n1-k1+1) + ...
          gammaln(n+1) - gammaln(k+1) - gammaln(n-k+1) + ...
          gammaln(k1+k+1) + gammaln(n+n1-k-k1+1) - gammaln(n+n1+2));
  p = max(0,min(1,p));
end