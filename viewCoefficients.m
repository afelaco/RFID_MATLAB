function viewCoefficients(c, t, max)

%% Degree and order vectors.
L = sqrt(length(c))-1;
i = 0:(L+1)^2-1;
l = floor(sqrt(i));
m = i-l.*(l+1);

c = cell2mat(c);

dim = size(c, 2);

for i = 0 : L
    for j = 1 : dim
        
        C(i+1, L-i+1:L+i+1, j) = c(l == i, j);
        
    end
end

l = 0:L;
m = -L:L;

figure
for i = 1 : dim
    
    clims = [0 max];
    
    subplot(dim, 1, i)
    imagesc(m,l,abs(C(:,:,i)),'AlphaData',l' >= abs(m),clims)
    title(t(i), 'interpreter', 'latex')
    axis equal
    xlim([-L-0.5 L+0.5])
    ylim([-0.5 L+0.5])
    colorbar
    
end

end