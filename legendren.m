function output = legendren(l, m, x)

output = zeros(size(x));

if abs(m) <= l
    
    output = legendre(l, x);
    
    output = output(abs(m)+1, :);
    
    if m < 0
        
        output = (-1)^m.*factorial(l+m)./factorial(l-m).*output;
        
    end
    
end

end