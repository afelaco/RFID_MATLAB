function output = ALF(l, x)

output = legendre(l, x);

if l > 0

    output = permute(output, [2 3 1]);
    
end

end