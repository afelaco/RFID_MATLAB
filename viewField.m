function viewField(Field, Grid)

[x, y, z] = sph2cart_c(Grid.Theta, Grid.Phi, Grid.Rho);

xmax = max(Grid.Rho, [], 'all')/sqrt(2);

xq = linspace(-xmax, xmax);

[xq, yq, zq] = meshgrid(xq);

if ndims(Field) >= 4
    
    v = vecnorm(Field, 2, 4);
    
else
    
    v = abs(Field);
    
end

vq = griddata(x, y, z, v, xq, yq, zq);

figure
h = slice(xq, yq, zq, vq, 0, 0, 0);
axis equal
set(h, 'edgecolor', 'none')
set(h, 'facecolor', 'interp')
view([1 1 1])

end