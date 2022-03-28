function [boundary_dof_back_face] = get_boundary_dof_back2(dtize_basis,dtize_nodes, Kx_coarse, Ky_coarse, Kz_coarse, Kx_fine, Ky_fine, Kz_fine, test_phi, dtize_weights, element_nodes2, domain_mapping)

%% back face

basis_back_face = kron(dtize_basis.ey,dtize_basis.ex);

[nodes_back_face.eta,nodes_back_face.xii,nodes_back_face.zta] = meshgrid(dtize_nodes.sx,dtize_nodes.sy,-1);

nodes_back_face.xii = nodes_back_face.xii(:)';
nodes_back_face.eta = nodes_back_face.eta(:)';
nodes_back_face.zta = nodes_back_face.zta(:)';

for i1 = 1:Kx_coarse
    for j1 = 1:Ky_coarse
        for k1 = 1
            eleid1 = (j1 - 1)*Kx_coarse + i1;
            for i2 = 1:Kx_fine
                for j2 = 1:Ky_fine
                    for k2 = 1
                        eleid2 = (eleid1 - 1)*Kx_fine*Ky_fine + (j2-1)*Kx_fine + i2;
                        [nodes_back_face.xii2(:,:,eleid2),nodes_back_face.eta2(:,:,eleid2),nodes_back_face.zta2(:,:,eleid2)] = element_nodes2(nodes_back_face.xii,nodes_back_face.eta,nodes_back_face.zta,i1,j1,k1,i2,j2,k2);
                    end
                end
            end
        end
    end
end

[nodes_back_face.xxx,nodes_back_face.yyy,nodes_back_face.zzz] = domain_mapping(nodes_back_face.xii2,nodes_back_face.eta2,nodes_back_face.zta2);

phi_back_face = test_phi(nodes_back_face.xxx,nodes_back_face.yyy,nodes_back_face.zzz);

w88_back_face = kron(dtize_weights.wy,dtize_weights.wx);

nn_back = -1;

Kx = Kx_coarse * Kx_fine;
Ky = Ky_coarse * Ky_fine;
% Kz = Kz_coarse * Kz_fine;

for eln = 1:Kx*Ky
    boundary_dof_back_face(:,eln) = nn_back * reduction(phi_back_face(:,:,eln), basis_back_face, w88_back_face, 1);
end

end

