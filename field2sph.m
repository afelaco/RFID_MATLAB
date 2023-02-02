function field2sph(Field, Theta, Phi)

Rho = ones(size(Theta));

X = Rho .* sin(Theta) .* cos(Phi);
Y = Rho .* sin(Theta) .* sin(Phi);
Z = Rho .* cos(Theta);

surf(X,Y,Z,Field,'EdgeColor','none')
axis equal

end