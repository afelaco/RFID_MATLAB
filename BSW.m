function J = BSW(k, Grid, Y)

%% Subgrid.
Rho = Grid.Rho(1,1,:);

%% Order.
L = sqrt(length(Y))-1;

for l = 0 : L
    for m = -l : l
   
        J(:, :, :, l*(l+1)+m+1) = (1i)^(-l).*sphbsl(l, k, Rho).*extract(Y, l, m);
        
    end
end

%% Pack HSW.
J = squeeze(num2cell(J, [1 2 3]));

end