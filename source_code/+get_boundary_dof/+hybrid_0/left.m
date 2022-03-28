function [boundary_dof_left_face] = get_hybrid_0_boundary_dof_left(Kx,Ky,Kz,p,dtize_basis,dtize_nodes,dtize_weights,element_nodes,domain_mapping,test_phi)

%% left face

basis_left_face = kron(dtize_basis.ez,dtize_basis.ey);

[nodes_left_face.xii,nodes_left_face.eta,nodes_left_face.zta] = meshgrid(-1,dtize_nodes.sy,dtize_nodes.sz);

nodes_left_face.xii = nodes_left_face.xii(:)';
nodes_left_face.eta = nodes_left_face.eta(:)';
nodes_left_face.zta = nodes_left_face.zta(:)';

% get new ref nodes
[nodes_left_face.xii2,nodes_left_face.eta2,nodes_left_face.zta2] = new_ref_nodes(1,1:Ky,1:Kz,nodes_left_face.xii,nodes_left_face.eta,nodes_left_face.zta,Ky,element_nodes);

[nodes_left_face.xxx,nodes_left_face.yyy,nodes_left_face.zzz] = domain_mapping(nodes_left_face.xii2,nodes_left_face.eta2,nodes_left_face.zta2);

phi_left_face = test_phi(nodes_left_face.xxx,nodes_left_face.yyy,nodes_left_face.zzz);

w88_left_face = kron(dtize_weights.wz,dtize_weights.wy);

nn_left = -1;

boundary_dof_left_face = zeros(p^2,Ky*Kz);

for eln = 1:Ky*Kz
    boundary_dof_left_face(:,eln) = nn_left * reduction(phi_left_face(:,:,eln), basis_left_face, w88_left_face, 1);
end

end

