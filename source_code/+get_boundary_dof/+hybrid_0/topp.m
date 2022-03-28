function [boundary_dof_topp_face] = get_hybrid_0_boundary_dof_topp(Kx,Ky,Kz,~,dtize_basis,dtize_nodes,dtize_weights,element_nodes,domain_mapping,test_phi)

%% top face

basis_topp_face = kron(dtize_basis.ez,dtize_basis.ex);

[nodes_topp_face.xii,nodes_topp_face.eta,nodes_topp_face.zta] = meshgrid(dtize_nodes.sx,+1,dtize_nodes.sz);

nodes_topp_face.xii = nodes_topp_face.xii(:)';
nodes_topp_face.eta = nodes_topp_face.eta(:)';
nodes_topp_face.zta = nodes_topp_face.zta(:)';

for i = 1:Kx
    for k = 1:Kz
        local_eleid = (k-1)*Kx + i;
        [dtize.nodes_topp_face.xii2(:,:,local_eleid),dtize.nodes_topp_face.eta2(:,:,local_eleid),dtize.nodes_topp_face.zta2(:,:,local_eleid)] = element_nodes(nodes_topp_face.xii,nodes_topp_face.eta,nodes_topp_face.zta,i,Ky,k);
    end
end

w88_topp_face = kron(dtize_weights.wz,dtize_weights.wx);

nn_topp = +1;

boundary_dof_topp_face = generate_boundary_dof(dtize.nodes_topp_face, domain_mapping, test_phi, basis_topp_face, w88_topp_face, nn_topp, Kx*Kz);

end

