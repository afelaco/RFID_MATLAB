function output = SRF(l, k, x)

Xj = besselzero(l+0.5, 1, 1);

Xy = besselzero(l+0.5, 1, 2);

output(x < Xy/k) = sphbsl(l, k.*x(x < Xy/k)); % + 1i*sphbsl(l, Xj./Xy.*k.*x(x < Xy/k));

output(x >= Xy/k) = sphhnk(l, 2, k.*x(x >= Xy/k));

output = reshape(output, size(x));

end