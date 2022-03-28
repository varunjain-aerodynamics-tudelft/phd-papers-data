function [boundary_dof_rght_face] = get_hybrid_0_boundary_dof_rght(Kx,Ky,Kz,p,dtize_basis,dtize_nodes,dtize_weights,element_nodes,domain_mapping,test_phi)

%% right face

basis_rght_face = kron(dtize_basis.ez,dtize_basis.ey);

[nodes_rght_face.xii,nodes_rght_face.eta,nodes_rght_face.zta] = meshgrid(+1,dtize_nodes.sy,dtize_nodes.sz);

nodes_rght_face.xii = nodes_rght_face.xii(:)';
nodes_rght_face.eta = nodes_rght_face.eta(:)';
nodes_rght_face.zta = nodes_rght_face.zta(:)';

% get new ref nodes
[nodes_rght_face.xii2,nodes_rght_face.eta2,nodes_rght_face.zta2] = new_ref_nodes(Kx,1:Ky,1:Kz,nodes_rght_face.xii,nodes_rght_face.eta,nodes_rght_face.zta,Ky,element_nodes);

[nodes_rght_face.xxx,nodes_rght_face.yyy,nodes_rght_face.zzz] = domain_mapping(nodes_rght_face.xii2,nodes_rght_face.eta2,nodes_rght_face.zta2);

phi_rght_face = test_phi(nodes_rght_face.xxx,nodes_rght_face.yyy,nodes_rght_face.zzz);

w88_rght_face = kron(dtize_weights.wz,dtize_weights.wy);

nn_rght = +1;

for eln = 1:Ky*Kz
    boundary_dof_rght_face(:,eln) = nn_rght * reduction(phi_rght_face(:,:,eln), basis_rght_face, w88_rght_face, 1); 
end

end

