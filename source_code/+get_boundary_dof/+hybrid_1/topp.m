function [boundary_dof_topp_face] = get_boundary_dof_topp2(dtize_basis,dtize_nodes, Kx_coarse, Ky_coarse, Kz_coarse, Kx_fine, Ky_fine, Kz_fine, test_phi, dtize_weights, element_nodes2, domain_mapping)

%% top face

basis_topp_face = kron(dtize_basis.ez,dtize_basis.ex);

[nodes_topp_face.xii,nodes_topp_face.eta,nodes_topp_face.zta] = meshgrid(dtize_nodes.sx,+1,dtize_nodes.sz);

nodes_topp_face.xii = nodes_topp_face.xii(:)';
nodes_topp_face.eta = nodes_topp_face.eta(:)';
nodes_topp_face.zta = nodes_topp_face.zta(:)';

for i1 = 1:Kx_coarse
    for j1 = Ky_coarse
        for k1 = 1:Kz_coarse
            eleid1 = (k1 - 1)*Kx_coarse + i1;
            for i2 = 1:Kx_fine
                for j2 = Ky_fine
                    for k2 = 1:Kz_fine
                        eleid2 = (eleid1 - 1)*Kx_fine*Kz_fine + (k2-1)*Kx_fine + i2;
                        [nodes_topp_face.xii2(:,:,eleid2),nodes_topp_face.eta2(:,:,eleid2),nodes_topp_face.zta2(:,:,eleid2)] = element_nodes2(nodes_topp_face.xii,nodes_topp_face.eta,nodes_topp_face.zta,i1,j1,k1,i2,j2,k2);
                    end
                end
            end
        end
    end
end

[nodes_topp_face.xxx,nodes_topp_face.yyy,nodes_topp_face.zzz] = domain_mapping(nodes_topp_face.xii2,nodes_topp_face.eta2,nodes_topp_face.zta2);

phi_topp_face = test_phi(nodes_topp_face.xxx,nodes_topp_face.yyy,nodes_topp_face.zzz);

w88_topp_face = kron(dtize_weights.wz,dtize_weights.wx);

nn_topp = +1;

Kx = Kx_coarse * Kx_fine;
% Ky = Ky_coarse * Ky_fine;
Kz = Kz_coarse * Kz_fine;

for eln = 1:Kx*Kz
    boundary_dof_topp_face(:,eln) = nn_topp * reduction(phi_topp_face(:,:,eln), basis_topp_face, w88_topp_face, 1);
end

end

