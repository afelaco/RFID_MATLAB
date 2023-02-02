function output = sphhnk(l, K, k, r)

x = k.*r;

output = sqrt(pi./(2.*x)).*besselh(l+0.5, K, x);

end