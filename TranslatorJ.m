function T = TranslatorJ(s, t, p, k0, L)

for la = 0 : L
    for mu = -la : la
        
        i = la*(la+1)+mu+1;
        
        for l = 0 : L
            for m = -l : l
                
                j = l*(l+1)+m+1;
                
                for q = 0 : L
                    if abs(m-mu) <= q
                        
                        k = q+1;
                        
                        T(i, j, k) = 4*pi*(-1)^(la-l+mu)*Gaunt(l, m, la, -mu, q).* ...
                            sqrt((2*q+1)./(4*pi).*factorial(q-m+mu)./factorial(q+m-mu)).* ...
                            (1i)^(-q).*sphbsl(q, k0*s)*legendren(q, m-mu, cos(t))*exp(-1i*(m-mu)*p);
                        
                    end
                end
            end
        end
    end
end

T = sum(T, 3);

end