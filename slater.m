function slater(L)

% Support Variables.
full_vec = (L+1)^2;

% Preallocation.
MM = zeros((L+1)^2,(L+1)^2,(L+1)^2);
ME = MM;
counter = 0;

% Degree and Order Vectors.
index = (1:full_vec)';
l = ceil(sqrt(index))-1;
m = index-l.*(l+1)-1;

% Slater Matrix.
for i = 1 : full_vec
    for j = 1 : full_vec
        for k = 1 : full_vec
            
            l_1 = l(i);
            l_2 = l(j);
            l_3 = l(k);
            
            m_1 = m(i);
            m_2 = m(j);
            m_3 = m(k);
            
            if m_1 + m_2 == m_3
                
                if rem(l_1+l_2+l_3, 2) == 0
                    
                    if l_3 >= abs(l_1-l_2) && l_3 <= l_1+l_2
                        
                        for m_s = max([-1, m_2-l_2, m_3-l_3]) : min([1, m_2+l_2, m_3+l_3])
                            
                            counter(end+1) = CG(l_2,m_2-m_s,1,m_s,l_2,m_2).*CG(l_3,m_3-m_s,1,m_s,l_3,m_3).*CG(l_1,m_1,l_2,m_2-m_s,l_3,m_3-m_s);
                            
                        end
                        
                        MM(i,j,k) = (-1).^l_3.*sqrt((2.*l_1+1).*(2.*l_2+1)./(2.*l_3+1)./(4.*pi)).*CG(l_1,0,l_2,0,l_3,0).*sum(counter);
                        
                        counter = 0;
                        
                    end
                    
                else
                    
                    if l_3 > abs(l_1-l_2) && l_3 <= l_1+l_2+1
                        
                        for m_s = max([-1, m_2-l_2, m_3-l_3+1]) : min([1, m_2 + l_2, m_3+l_3-1])
                            
                            counter(end+1) = CG(l_2,m_2-m_s,1,m_s,l_2,m_2).*CG(l_3-1,m_3-m_s,1,m_s,l_3,m_3).*CG(l_1,m_1,l_2,m_2-m_s,l_3-1,m_3-m_s);
                            
                        end
                        
                        ME(i,j,k) = sqrt((l_3+1)./(2.*l_3+1)).*(-1).^l_3.*sqrt((2.*l_1+1).*(2.*l_2+1)./(2.*l_3-1)./(4.*pi)).*CG(l_1,0,l_2,0,l_3-1,0).*sum(counter);
                        
                        counter = 0;
                        
                    end
                    
                    if l_3 >= abs(l_1-l_2)-1 && l_3 < l_1+l_2
                        
                        for m_s = max([-1, m_2-l_2, m_3-l_3-1]) : min([1, m_2+l_2, m_3+l_3+1])
                            
                            counter(end+1) = CG(l_2,m_2-m_s,1,m_s,l_2,m_2).*CG(l_3+1,m_3-m_s,1,m_s,l_3,m_3).*CG(l_1,m_1,l_2,m_2-m_s,l_3+1,m_3-m_s);
                            
                        end
                      
                        ME(i,j,k) = ME(i,j,k) + sqrt((l_3)./(2.*l_3+1)).*(-1).^l_3.*sqrt((2.*l_1+1).*(2.*l_2+1)./(2.*l_3+3)./(4.*pi)).*CG(l_1,0,l_2,0,l_3+1,0).*sum(counter);
                        
                        counter = 0;
                        
                    end
                end
            end
        end
    end
end

save('Lib\Slater','MM','ME')

end