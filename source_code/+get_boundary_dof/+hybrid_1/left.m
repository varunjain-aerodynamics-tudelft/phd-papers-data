function [boundary_dof_left_face] = get_boundary_dof_left2(dtize_basis,dtize_nodes, Kx_coarse, Ky_coarse, Kz_coarse, Kx_fine, Ky_fine, Kz_fine, test_phi, dtize_weights, element_nodes2, domain_mapping)

%% left face

basis_left_face = kron(dtize_basis.ez,dtize_basis.ey);

[nodes_left_face.xii,nodes_left_face.eta,nodes_left_face.zta] = meshgrid(-1,dtize_nodes.sy,dtize_nodes.sz);

nodes_left_face.xii = nodes_left_face.xii(:)';
nodes_left_face.eta = nodes_left_face.eta(:)';
nodes_left_face.zta = nodes_left_face.zta(:)';

% get new ref nodes
% [nodes_left_face.xii2,nodes_left_face.eta2,nodes_left_face.zta2] = new_ref_nodes(1,1:Ky,1:Kz,nodes_left_face.xii,nodes_left_face.eta,nodes_left_face.zta,Ky,element.nodes);

% [nodes_left_face.xxx,nodes_left_face.yyy,nodes_left_face.zzz] = domain.mapping(nodes_left_face.xii2,nodes_left_face.eta2,nodes_left_face.zta2);
% ttl_nr_sub_el = Kx_fine*Ky_fine*Kz_fine;

% for eln1 = 1:Ky_coarse*Kz_coarse
%     for eln2 = 1:Ky_fine*Kz_fine
%         eln = (mac_ele_id_left(eln1)-1)*ttl_nr_sub_el + sub_ele_id_left(eln2);
%         eleid = (eln1-1)*Ky_fine*Kz_fine + eln2;
%         [nodes_left_face.xxx(:,:,eleid),nodes_left_face.yyy(:,:,eleid),nodes_left_face.zzz(:,:,eleid)] = mesh_fine.el(eln).mapping(nodes_left_face.xii,nodes_left_face.eta,nodes_left_face.zta);
%     end
% end

for i1 = 1
    for j1 = 1:Ky_coarse
        for k1 = 1:Kz_coarse
            eleid1 = (k1 - 1)*Ky_coarse + j1;
%             eleid_coarse = (k1-1)*Kx_coarse*Ky_coarse + (j1-1)*Kx_coarse + i1;
            for i2 = 1
                for j2 = 1:Ky_fine
                    for k2 = 1:Kz_fine
%                         eleid_fine = (eleid_coarse-1)*ttl_nr_sub_el + (k2-1)*Kx_fine*Ky_fine + (j2-1)*Kx_fine + i2;
                        eleid2 = (eleid1 - 1)*Ky_fine*Kz_fine + (k2-1)*Ky_fine + j2;
                        [nodes_left_face.xii2(:,:,eleid2),nodes_left_face.eta2(:,:,eleid2),nodes_left_face.zta2(:,:,eleid2)] = element_nodes2(nodes_left_face.xii,nodes_left_face.eta,nodes_left_face.zta,i1,j1,k1,i2,j2,k2);
                    end
                end
            end
        end
    end
end

[nodes_left_face.xxx,nodes_left_face.yyy,nodes_left_face.zzz] = domain_mapping(nodes_left_face.xii2,nodes_left_face.eta2,nodes_left_face.zta2);

phi_left_face = test_phi(nodes_left_face.xxx,nodes_left_face.yyy,nodes_left_face.zzz);

w88_left_face = kron(dtize_weights.wz,dtize_weights.wy);

nn_left = -1;

%%%%%%% ERRORRRR here !!!!! 

% Kx = Kx_coarse * Kx_fine;
Ky = Ky_coarse * Ky_fine;
Kz = Kz_coarse * Kz_fine;

for eln = 1:Ky*Kz
    boundary_dof_left_face(:,eln) = nn_left * reduction(phi_left_face(:,:,eln), basis_left_face, w88_left_face, 1);
end

end

