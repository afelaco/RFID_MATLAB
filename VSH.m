function [E, M] = VSH(Y)

%% Order.
L = sqrt(length(Y))-1;

for l = 0 : L
    for m = -l : l
        for ms = -1 : 1
            
            E(:, :, mod(ms+2,3)+1, l*(l+1)+m+1) = sqrt((l+1)/(2*l+1)).*CG(l-1, m-ms, 1, ms, l, m).*extract(Y, l-1, m-ms) + ...
                sqrt(l/(2*l+1)).*CG(l+1, m-ms, 1, ms, l, m).*extract(Y, l+1, m-ms);
            
            M(:, :, mod(ms+2,3)+1, l*(l+1)+m+1) = CG(l, m-ms, 1, ms, l, m).*extract(Y, l, m-ms);
            
        end
    end
end

%% Pack VSH.
E = squeeze(num2cell(E, [1 2 3]));
M = squeeze(num2cell(M, [1 2 3]));

end