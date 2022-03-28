function [mac_ele_id] = get_boundary_element_ID(Kx_coarse,Ky_coarse, Kz_coarse)

%%%%%%%%%%%%%%%% left %%%%%%%%%%%%%%%%%

for j = 1:Ky_coarse
    for k = 1:Kz_coarse
        local_id = (k-1)*Ky_coarse + j;
        globl_id = (k-1)*Kx_coarse*Ky_coarse + (j-1)*Kx_coarse + 1;
        mac_ele_id.left(local_id) = globl_id;
    end
end

%%%%%%%%%%%%%%%% right %%%%%%%%%%%%%%%%%

for j = 1:Ky_coarse
    for k = 1:Kz_coarse
        local_id = (k-1)*Ky_coarse + j;
        globl_id = (k-1)*Kx_coarse*Ky_coarse + (j-1)*Kx_coarse + Kx_coarse;
        mac_ele_id.rght(local_id) = globl_id;
    end
end

%%%%%%%%%%%%%%%% bottom %%%%%%%%%%%%%%%%%

for i = 1:Kx_coarse
    for k = 1:Kz_coarse
        local_id = (k-1)*Kx_coarse + i;
        globl_id = (k-1)*Kx_coarse*Ky_coarse + (1-1)*Kx_coarse + i;
        mac_ele_id.bttm(local_id) = globl_id;
    end
end

%%%%%%%%%%%%%%%% top %%%%%%%%%%%%%%%%%

for i = 1:Kx_coarse
    for k = 1:Kz_coarse
        local_id = (k-1)*Kx_coarse + i;
        globl_id = (k-1)*Kx_coarse*Ky_coarse + (Ky_coarse-1)*Kx_coarse + i;
        mac_ele_id.topp(local_id) = globl_id;
    end
end

%%%%%%%%%%%%%%%% back %%%%%%%%%%%%%%%%%

for i = 1:Kx_coarse
    for j = 1:Ky_coarse
        local_id = (j-1)*Kx_coarse + i;
        globl_id = (1-1)*Kx_coarse*Ky_coarse + (j-1)*Kx_coarse + i;
        mac_ele_id.back(local_id) = globl_id;
    end
end

%%%%%%%%%%%%%%%% front %%%%%%%%%%%%%%%%%

for i = 1:Kx_coarse
    for j = 1:Ky_coarse
        local_id = (j-1)*Kx_coarse + i;
        globl_id = (Kz_coarse-1)*Kx_coarse*Ky_coarse + (j-1)*Kx_coarse + i;
        mac_ele_id.frnt(local_id) = globl_id;
    end
end

end

