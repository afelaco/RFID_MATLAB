function R = RSW(k, Grid, Y)

%% Subgrid.
Rho = Grid.Rho(1,1,:);

%% Order.
L = sqrt(length(Y))-1;

for l = 0 : L
    for m = -l : l
   
        R(:, :, :, l*(l+1)+m+1) = (1i)^(-l).*sphrdl(l, k, Rho).*extract(Y, l, m);
        
    end
end

%% Pack HSW.
R = squeeze(num2cell(R, [1 2 3]));

end