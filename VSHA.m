function [e, m] = VSHA(FarField, L)

%% Subgrid.
[Grid, order] = GridLebedev(2*L);

%% SSH.
Tag = strcat('Lib\SSH\', sprintf('lebedev_%d.mat', order));

if isfile(Tag) == 0
    
    Y = SSH(L, Grid);
    save(Tag, 'Y')
    
else
    
    load(Tag)
    
end

%% VSH.
Tag = strcat('Lib\VSH\', sprintf('lebedev_%d.mat', order));

if isfile(Tag) == 0
    
    [E, M] = VSH(Y);
    save(Tag, 'E', 'M')
    
else
    
    load(Tag)
    
end

%% Interpolate data.
FieldPol(:,:,1) = interp2(FarField.Phi, FarField.Theta, FarField.Field(:,:,1), Grid.Phi, Grid.Theta, 'spline');
FieldPol(:,:,2) = interp2(FarField.Phi, FarField.Theta, FarField.Field(:,:,2), Grid.Phi, Grid.Theta, 'spline');

%% Change of basis.
FieldSph(:,:,1) = exp(-1i.*Grid.Phi)./sqrt(2).*(-cos(Grid.Theta).*FieldPol(:,:,1) + 1i.*FieldPol(:,:,2));
FieldSph(:,:,2) = exp(1i.*Grid.Phi)./sqrt(2).*(cos(Grid.Theta).*FieldPol(:,:,1) + 1i.*FieldPol(:,:,2));
FieldSph(:,:,3) = -sin(Grid.Theta).*FieldPol(:,:,1);

%% Convolution.
e = squeeze(sum(Grid.Weight.*FieldSph.*conj(cat(4, E{:})), [1 2 3]));
m = squeeze(sum(Grid.Weight.*FieldSph.*conj(cat(4, M{:})), [1 2 3]));

end