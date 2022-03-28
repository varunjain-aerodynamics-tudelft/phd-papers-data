function [boundary_global_ID] = get_hybrid_1_boundary_element_ID(Kx_coarse,Ky_coarse,Kz_coarse,Kx_fine,Ky_fine,Kz_fine)

mac_bndry_ele_id = mesh.dim_3.get_boundary_element_ID(Kx_coarse,Ky_coarse, Kz_coarse);
sub_bndry_ele_id = mesh.dim_3.get_boundary_element_ID(Kx_fine,Ky_fine, Kz_fine);

ttl_nr_sub_el = Kx_fine * Ky_fine * Kz_fine;

%%%%%%%%%%%%%%%% left %%%%%%%%%%%%%%%%%

for eln1 = 1:Ky_coarse*Kz_coarse
    for eln2 = 1:Ky_fine*Kz_fine
        local_eleid = (eln1-1)*Ky_fine*Kz_fine + eln2;
        boundary_global_ID.left(local_eleid) = (mac_bndry_ele_id.left(eln1)-1)*ttl_nr_sub_el + sub_bndry_ele_id.left(eln2);
    end
end

%%%%%%%%%%%%%%%% right %%%%%%%%%%%%%%%%%

for eln1 = 1:Ky_coarse*Kz_coarse
    for eln2 = 1:Ky_fine*Kz_fine
        local_eleid = (eln1-1)*Ky_fine*Kz_fine + eln2;
        boundary_global_ID.rght(local_eleid) = (mac_bndry_ele_id.rght(eln1)-1)*ttl_nr_sub_el + sub_bndry_ele_id.rght(eln2);
    end
end

%%%%%%%%%%%%%%%% bottom %%%%%%%%%%%%%%%%%

for eln1 = 1:Kx_coarse*Kz_coarse
    for eln2 = 1:Kx_fine*Kz_fine
        local_eleid = (eln1-1)*Kx_fine*Kz_fine + eln2;
        boundary_global_ID.bttm(local_eleid) = (mac_bndry_ele_id.bttm(eln1)-1)*ttl_nr_sub_el + sub_bndry_ele_id.bttm(eln2);
    end
end

%%%%%%%%%%%%%%%% top %%%%%%%%%%%%%%%%%

for eln1 = 1:Kx_coarse*Kz_coarse
    for eln2 = 1:Kx_fine*Kz_fine
        local_eleid = (eln1-1)*Kx_fine*Kz_fine + eln2;
        boundary_global_ID.topp(local_eleid) = (mac_bndry_ele_id.topp(eln1)-1)*ttl_nr_sub_el + sub_bndry_ele_id.topp(eln2);
    end
end

%%%%%%%%%%%%%%%% back %%%%%%%%%%%%%%%%%

for eln1 = 1:Kx_coarse*Ky_coarse
    for eln2 = 1:Kx_fine*Ky_fine
        local_eleid = (eln1-1)*Kx_fine*Ky_fine + eln2;
        boundary_global_ID.back(local_eleid) = (mac_bndry_ele_id.back(eln1)-1)*ttl_nr_sub_el + sub_bndry_ele_id.back(eln2);
    end
end

%%%%%%%%%%%%%%%% front %%%%%%%%%%%%%%%%%

for eln1 = 1:Kx_coarse*Ky_coarse
    for eln2 = 1:Kx_fine*Ky_fine
        local_eleid = (eln1-1)*Kx_fine*Ky_fine + eln2;
        boundary_global_ID.frnt(local_eleid) = (mac_bndry_ele_id.frnt(eln1)-1)*ttl_nr_sub_el + sub_bndry_ele_id.frnt(eln2);
    end
end

end

