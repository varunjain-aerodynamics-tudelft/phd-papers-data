function [boundary_dof_bttm_face] = generate_boundary_dof(dtize_nodes_bttm_face, domain_mapping, test_phi, basis_bttm_face, w88_bttm_face, nn_bttm, nr_face_el)

% basis_bttm_face = kron(dtize.basis.ez,dtize.basis.ex);
% [nodes_bttm_face.xii,nodes_bttm_face.eta,nodes_bttm_face.zta] = meshgrid(dtize.nodes.sx,-1,dtize.nodes.sz);

% nodes_bttm_face.xii = nodes_bttm_face.xii(:)';
% nodes_bttm_face.eta = nodes_bttm_face.eta(:)';
% nodes_bttm_face.zta = nodes_bttm_face.zta(:)';
% 
% % get new ref nodes
% 
% for i = 1:Kx
%     for k = 1:Kz
%         local_eleid = (k-1)*Kx + i;
%         [dtize.nodes_bttm_face.xii2(:,:,local_eleid),dtize.nodes_bttm_face.eta2(:,:,local_eleid),dtize.nodes_bttm_face.zta2(:,:,local_eleid)] = element.nodes(nodes_bttm_face.xii,nodes_bttm_face.eta,nodes_bttm_face.zta,i,1,k);
%     end
% end

[nodes_bttm_face.xxx,nodes_bttm_face.yyy,nodes_bttm_face.zzz] = domain_mapping(dtize_nodes_bttm_face.xii2,dtize_nodes_bttm_face.eta2,dtize_nodes_bttm_face.zta2);

phi_bttm_face = test_phi(nodes_bttm_face.xxx,nodes_bttm_face.yyy,nodes_bttm_face.zzz);

% w88_bttm_face = kron(dtize.weights.wz,dtize.weights.wx);

% nn_bttm = -1;

for eln = 1:nr_face_el
    boundary_dof_bttm_face(:,eln) = nn_bttm * reduction(phi_bttm_face(:,:,eln), basis_bttm_face, w88_bttm_face, 1);
end

end

