function hM = hMatrix(L)

i = (0 : (L+1)^2-1)';
l = floor(sqrt(i));
m = i - l.*(l+1);
ms = reshape([+1, -1, 0], 1, 1, []);

hM.e = (l+1 == l').*(m+ms == m').*(abs(m+ms) <= l').*sqrt((l+2)./(2.*l+3)) + ...
    (l-1 == l').*(m+ms == m').*(abs(m+ms) <= l').*sqrt((l-1)./(2.*l-1));

hM.m = (l == l').*(m+ms == m').*(abs(m+ms) <= l');

for l1 = 0 : L
    for m1 = -l1 : l1
        for l3 = 0 : L
            for m3 = -l3 : l3
                for ms = -1 : 1
                    
                    CGM(l1*(l1+1)+m1+1, l3*(l3+1)+m3+1, mod(ms+2,3)+1) = CG(l1, m1, 1, ms, l3, m3);
                    
                end
            end
        end
    end
end

hM.e = CGM.*hM.e;
hM.m = CGM.*hM.m;
hM.L = L;

end