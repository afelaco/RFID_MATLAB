function Y = SSH(L, Grid)

%% Subgrid.
Theta = Grid.Theta(:,:,1);
Phi = Grid.Phi(:,:,1);

%% Preallocation.
Y = zeros([size(Theta), (L+1)^2]);

%% Degree and order vectors.
i = reshape(0:(L+1)^2-1, 1, 1, []);
l = floor(sqrt(i));
m = i-l.*(l+1);

%% Associated Legendre functions.
for i = 0 : L
    
    Y(:,:,l==i & m>=0) = ALF(i,cos(Theta));
    
end

%% Phase.
Y = Y.*exp(1i.*m.*Phi);

%% Normalization.
Y = Y.*sqrt((2.*l+1)./(4*pi).*factorial(l-m)./factorial(l+m));

%% Recursion.
for i = 0 : L
   
    Y(:,:,l == i & m < 0) = flip((-1).^m(:,:,l == i & m > 0).*conj(Y(:,:,l == i & m > 0)), 3);
    
end

%% Pack.
Y = squeeze(num2cell(Y, [1 2]));

end