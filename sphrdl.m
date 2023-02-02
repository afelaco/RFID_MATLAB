function output = sphrdl(l, k, r)

output = sphbsl(l, k, r) + 1i.*sphbsl(l-1, k, r);

end