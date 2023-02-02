function H = HSW(k, Grid, Y)

%% Subgrid.
Rho = Grid.Rho(1,1,:);

%% Order.
L = sqrt(length(Y))-1;

for l = 0 : L
    for m = -l : l
   
        H(:, :, :, l*(l+1)+m+1) = (1i)^(-l).*sphhnk(l, 2, k, Rho).*extract(Y, l, m);
        
    end
end

%% Pack HSW.
H = squeeze(num2cell(H, [1 2 3]));

end