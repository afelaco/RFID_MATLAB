function output = Gaunt(l, m, la, mu, q)

output = (-1)^(m+mu)*sqrt((2*l+1)*(2*la+1)*(2*q+1)/(4*pi))* ...
    Wigner3j(l, 0, la, 0, q, 0)*Wigner3j(l, m, la, mu, q, -m-mu);

end