function [E10] = primal_grad_E10(p)

for i = 1:p
    for j = 1:p+1
        edgeij = (j-1)*p + i;
        
        node_left = (j-1)*(p+1) + i;
        node_rght = (j-1)*(p+1) + i+1;
        
        E10(edgeij,node_left) = +1;
        E10(edgeij,node_rght) = -1;
    end
end

for i = 1:p+1
    for j = 1:p
        edgeij = p*(p+1) + (j-1)*(p+1) + i;
        
        node_up = j*(p+1) + i;
        node_dn = (j-1)*(p+1) + i;
        
        E10(edgeij,node_up) = +1;
        E10(edgeij,node_dn) = -1;
    end
end

end

