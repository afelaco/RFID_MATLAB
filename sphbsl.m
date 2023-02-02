function output = sphbsl(l, k, r)

x = k.*r;

output = sqrt(pi./(2.*x)).*besselj(l+0.5, x);

if l == 0
    
    output(x == 0) = 1;
    
else
    
    output(x == 0) = 0;
    
end

output = reshape(output, size(x));

end