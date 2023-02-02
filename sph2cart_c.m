function [x, y, z] = sph2cart_c(r, t, p)

x = r .* cos(p) .* sin(t);
y = r .* sin(p) .* sin(t);
z = r .* cos(t);

end