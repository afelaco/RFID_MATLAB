function output = extract(input, l, m)

% Degree and order vectors.
L = sqrt(length(input))-1;

if l <= L && abs(m) <= l
    
    output = input{l*(l+1)+m+1};
    
else
    
    output = zeros(size(input{1}));
    
end

end