function D = WignerD(a, b, g, L)

%% Wigner-d matrix.
for l = 0 : L
    for m = -l : l
        for n = -l : l
   
            i = l*(l+1)+n+1;
            j = l*(l+1)+m+1;
            
            s = max(0, m-n):min(l+m, l-n);
            
            d(i,j) = sum((-1).^(n-m+s).*(cos(b./2)).^(2*(l-s)+m-n).*(sin(b./2)).^(n-m+2.*s)./ ...
                (factorial(l+m-s).*factorial(s).*factorial(n-m+s).*factorial(l-n-s)), 2);
            
        end
    end
end

%% Degree and order vectors.
i = (0:(L+1)^2-1)';
l = floor(sqrt(i));
m = i-l.*(l+1);

%% Normalization.
N = sqrt(factorial(l+m).*factorial(l-m).*factorial(l+m)'.*factorial(l-m)');

d = N.*d;

%% Wigner-D matrix.
D = exp(1i.*m.*a).*d.*exp(1i.*m'.*g);

end