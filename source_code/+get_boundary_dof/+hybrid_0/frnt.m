function [boundary_dof_frnt_face] = get_hybrid_0_boundary_dof_frnt(Kx,Ky,Kz,~,dtize_basis,dtize_nodes,dtize_weights,element_nodes,domain_mapping,test_phi)

%% front face

basis_frnt_face = kron(dtize_basis.ey,dtize_basis.ex);

[nodes_frnt_face.eta,nodes_frnt_face.xii,nodes_frnt_face.zta] = meshgrid(dtize_nodes.sx,dtize_nodes.sy,+1);

nodes_frnt_face.xii = nodes_frnt_face.xii(:)';
nodes_frnt_face.eta = nodes_frnt_face.eta(:)';
nodes_frnt_face.zta = nodes_frnt_face.zta(:)';

for i = 1:Kx
    for j = 1:Ky
        local_eleid = (j-1)*Kx + i;
        [nodes_frnt_face.xii2(:,:,local_eleid),nodes_frnt_face.eta2(:,:,local_eleid),nodes_frnt_face.zta2(:,:,local_eleid)] = element_nodes(nodes_frnt_face.xii,nodes_frnt_face.eta,nodes_frnt_face.zta,i,j,Kz);
    end
end

w88_frnt_face = kron(dtize_weights.wy,dtize_weights.wx);

nn_frnt = +1;

boundary_dof_frnt_face = generate_boundary_dof(nodes_frnt_face, domain_mapping, test_phi, basis_frnt_face, w88_frnt_face, nn_frnt, Kx*Ky);

end

