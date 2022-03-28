function [boundary_dof_back_face] = get_hybrid_0_boundary_dof_back(Kx,Ky,~,~,dtize_basis,dtize_nodes,dtize_weights,element_nodes,domain_mapping,test_phi)

%% back face

basis_back_face = kron(dtize_basis.ey,dtize_basis.ex);

[nodes_back_face.eta,nodes_back_face.xii,nodes_back_face.zta] = meshgrid(dtize_nodes.sx,dtize_nodes.sy,-1);

nodes_back_face.xii = nodes_back_face.xii(:)';
nodes_back_face.eta = nodes_back_face.eta(:)';
nodes_back_face.zta = nodes_back_face.zta(:)';

for i = 1:Kx
    for j = 1:Ky
        local_eleid = (j-1)*Kx + i;
        [nodes_back_face.xii2(:,:,local_eleid),nodes_back_face.eta2(:,:,local_eleid),nodes_back_face.zta2(:,:,local_eleid)] = element_nodes(nodes_back_face.xii,nodes_back_face.eta,nodes_back_face.zta,i,j,1);
    end
end

w88_back_face = kron(dtize_weights.wy,dtize_weights.wx);

nn_back = -1;

boundary_dof_back_face = generate_boundary_dof(nodes_back_face, domain_mapping, test_phi, basis_back_face, w88_back_face, nn_back, Kx*Ky);

end

