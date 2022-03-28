function [boundary_dof_bttm_face] = get_hybrid_0_boundary_dof_bttm(Kx,~,Kz,~,dtize_basis,dtize_nodes,dtize_weights,element_nodes,domain_mapping,test_phi)

%% bottom face

basis_bttm_face = kron(dtize_basis.ez,dtize_basis.ex);
[nodes_bttm_face.xii,nodes_bttm_face.eta,nodes_bttm_face.zta] = meshgrid(dtize_nodes.sx,-1,dtize_nodes.sz);

nodes_bttm_face.xii = nodes_bttm_face.xii(:)';
nodes_bttm_face.eta = nodes_bttm_face.eta(:)';
nodes_bttm_face.zta = nodes_bttm_face.zta(:)';

for i = 1:Kx
    for k = 1:Kz
        local_eleid = (k-1)*Kx + i;
        [dtize.nodes_bttm_face.xii2(:,:,local_eleid),dtize.nodes_bttm_face.eta2(:,:,local_eleid),dtize.nodes_bttm_face.zta2(:,:,local_eleid)] = element_nodes(nodes_bttm_face.xii,nodes_bttm_face.eta,nodes_bttm_face.zta,i,1,k);
    end
end

w88_bttm_face = kron(dtize_weights.wz,dtize_weights.wx);

nn_bttm = -1;

boundary_dof_bttm_face = generate_boundary_dof(dtize.nodes_bttm_face, domain_mapping, test_phi, basis_bttm_face, w88_bttm_face, nn_bttm, Kx*Kz);

end

