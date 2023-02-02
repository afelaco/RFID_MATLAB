function [h, L] = SSWA(e, m, ratio)

%% Order.
L = sqrt(length(e))-1;

Tag = 'Lib\hM';

if isfile(Tag) == 0
    
    hM = hMatrix(L);
    save(Tag, 'hM')
    
else
    
    load(Tag)
    
    if hM.L < L
        
        hM = hMatrix(L);
        save(Tag, 'hM')
        
    elseif hM.L > L
        
        hM.e = hM.e(1:(L+1)^2, 1:(L+1)^2, :);
        hM.m = hM.m(1:(L+1)^2, 1:(L+1)^2, :);
        
    end
    
end

hSph = squeeze(sum(hM.e.*cat(2, e{:}) + hM.m.*cat(2, m{:}), 2));

%% Trim.
norm = vecnorm(hSph, 2, 2);
norm = norm./max(norm, [], 'all');

boolean = norm >= ratio;

L = floor(sqrt(nnz(boolean))-1);

hSph = hSph(1:(L+1)^2, :);

%% Basis.
hCart(:,1) = 1/sqrt(2).*(hSph(:,2) - hSph(:,1));
hCart(:,2) = -1i/sqrt(2).*(hSph(:,2) + hSph(:,1));
hCart(:,3) = hSph(:,3);

%% Pack.
h = num2cell(hCart, 2);

end