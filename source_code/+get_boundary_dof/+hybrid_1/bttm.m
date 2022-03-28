function [boundary_dof_bttm_face] = get_boundary_dof_bttm2(dtize_basis,dtize_nodes, Kx_coarse, Ky_coarse, Kz_coarse, Kx_fine, Ky_fine, Kz_fine, test_phi, dtize_weights, element_nodes2, domain_mapping)

%% bottom face

basis_bttm_face = kron(dtize_basis.ez,dtize_basis.ex);
[nodes_bttm_face.xii,nodes_bttm_face.eta,nodes_bttm_face.zta] = meshgrid(dtize_nodes.sx,-1,dtize_nodes.sz);

nodes_bttm_face.xii = nodes_bttm_face.xii(:)';
nodes_bttm_face.eta = nodes_bttm_face.eta(:)';
nodes_bttm_face.zta = nodes_bttm_face.zta(:)';

for i1 = 1:Kx_coarse
    for j1 = 1
        for k1 = 1:Kz_coarse
            eleid1 = (k1 - 1)*Kx_coarse + i1;
            for i2 = 1:Kx_fine
                for j2 = 1
                    for k2 = 1:Kz_fine
                        eleid2 = (eleid1 - 1)*Kx_fine*Kz_fine + (k2-1)*Kx_fine + i2;
                        [nodes_bttm_face.xii2(:,:,eleid2),nodes_bttm_face.eta2(:,:,eleid2),nodes_bttm_face.zta2(:,:,eleid2)] = element_nodes2(nodes_bttm_face.xii,nodes_bttm_face.eta,nodes_bttm_face.zta,i1,j1,k1,i2,j2,k2);
                    end
                end
            end
        end
    end
end

[nodes_bttm_face.xxx,nodes_bttm_face.yyy,nodes_bttm_face.zzz] = domain_mapping(nodes_bttm_face.xii2,nodes_bttm_face.eta2,nodes_bttm_face.zta2);

phi_bttm_face = test_phi(nodes_bttm_face.xxx,nodes_bttm_face.yyy,nodes_bttm_face.zzz);

w88_bttm_face = kron(dtize_weights.wz,dtize_weights.wx);

nn_bttm = -1;

Kx = Kx_coarse * Kx_fine;
% Ky = Ky_coarse * Ky_fine;
Kz = Kz_coarse * Kz_fine;

for eln = 1:Kx*Kz
    boundary_dof_bttm_face(:,eln) = nn_bttm * reduction(phi_bttm_face(:,:,eln), basis_bttm_face, w88_bttm_face, 1);
end

end

