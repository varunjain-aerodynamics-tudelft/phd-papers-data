function [boundary_dof_frnt_face] = get_boundary_dof_frnt2(dtize_basis,dtize_nodes, Kx_coarse, Ky_coarse, Kz_coarse, Kx_fine, Ky_fine, Kz_fine, test_phi, dtize_weights, element_nodes2, domain_mapping)

%% front face

basis_frnt_face = kron(dtize_basis.ey,dtize_basis.ex);

[nodes_frnt_face.eta,nodes_frnt_face.xii,nodes_frnt_face.zta] = meshgrid(dtize_nodes.sx,dtize_nodes.sy,+1);

nodes_frnt_face.xii = nodes_frnt_face.xii(:)';
nodes_frnt_face.eta = nodes_frnt_face.eta(:)';
nodes_frnt_face.zta = nodes_frnt_face.zta(:)';

for i1 = 1:Kx_coarse
    for j1 = 1:Ky_coarse
        for k1 = Kz_coarse
            eleid1 = (j1 - 1)*Kx_coarse + i1;
            for i2 = 1:Kx_fine
                for j2 = 1:Ky_fine
                    for k2 = Kz_fine
                        eleid2 = (eleid1 - 1)*Kx_fine*Ky_fine + (j2-1)*Kx_fine + i2;
                        [nodes_frnt_face.xii2(:,:,eleid2),nodes_frnt_face.eta2(:,:,eleid2),nodes_frnt_face.zta2(:,:,eleid2)] = element_nodes2(nodes_frnt_face.xii,nodes_frnt_face.eta,nodes_frnt_face.zta,i1,j1,k1,i2,j2,k2);
                    end
                end
            end
        end
    end
end

[nodes_frnt_face.xxx,nodes_frnt_face.yyy,nodes_frnt_face.zzz] = domain_mapping(nodes_frnt_face.xii2,nodes_frnt_face.eta2,nodes_frnt_face.zta2);

phi_frnt_face = test_phi(nodes_frnt_face.xxx,nodes_frnt_face.yyy,nodes_frnt_face.zzz);

w88_frnt_face = kron(dtize_weights.wy,dtize_weights.wx);

nn_frnt = +1;

Kx = Kx_coarse * Kx_fine;
Ky = Ky_coarse * Ky_fine;

for eln = 1:Kx*Ky
    boundary_dof_frnt_face(:,eln) = nn_frnt * reduction(phi_frnt_face(:,:,eln), basis_frnt_face, w88_frnt_face, 1);
end

end

