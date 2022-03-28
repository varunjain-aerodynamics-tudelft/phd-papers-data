function [N_parallel] = unit_tangential(p)

% ---------- bottom boundary 

for i = 1:p
    bbtm_bndry_id = i;
    topp_bndry_id =  i;
        
    bttm_node_id = i;    
    topp_node_id = p*(p+1) + i
end

for j = 1:p
    
    rght_node_id = (p+1)*i;
    left_node_id = ()
end


