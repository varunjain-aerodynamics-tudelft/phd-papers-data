function [dof_f] = get_dof_3D_volumes(test_darcy_fff,basis_volume,dtize_metric, dtize_nodes, nr_volum, ttl_nr_el)

M333 = zeros(nr_volum, nr_volum, ttl_nr_el);

for eln = 1:ttl_nr_el
    M333(:,:,eln) = construct_mass_matrix(basis_volume,basis_volume,dtize_metric.www,dtize_metric.invg(:,:,eln));
end

ff_int = test_darcy_fff(dtize_nodes.x,dtize_nodes.y,dtize_nodes.z);

dof_f = zeros(nr_volum,1,ttl_nr_el);

for eln = 1:ttl_nr_el
    dof_f(:,1,eln) = M333(:,:,eln) \ reduction(ff_int(:,:,eln), basis_volume, dtize_metric.www, 1);
end

end

