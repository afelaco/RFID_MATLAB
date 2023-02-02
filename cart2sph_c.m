function [r, t, p] = cart2sph_c(x, y, z)

r = sqrt(x.^2 + y.^2 + z.^2);
t = acos(z./r);
p = atan(y./x);

end