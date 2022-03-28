function [boundary_dof_rght_face] = get_boundary_dof_rght2(dtize_basis,dtize_nodes, Kx_coarse, Ky_coarse, Kz_coarse, Kx_fine, Ky_fine, Kz_fine, test_phi, dtize_weights, element_nodes2, domain_mapping)

%% right face

basis_rght_face = kron(dtize_basis.ez,dtize_basis.ey);

[nodes_rght_face.xii,nodes_rght_face.eta,nodes_rght_face.zta] = meshgrid(+1,dtize_nodes.sy,dtize_nodes.sz);

nodes_rght_face.xii = nodes_rght_face.xii(:)';
nodes_rght_face.eta = nodes_rght_face.eta(:)';
nodes_rght_face.zta = nodes_rght_face.zta(:)';

for i1 = Kx_coarse
    for j1 = 1:Ky_coarse
        for k1 = 1:Kz_coarse
            eleid1 = (k1 - 1)*Ky_coarse + j1;
            for i2 = Kx_fine
                for j2 = 1:Ky_fine
                    for k2 = 1:Kz_fine
                        eleid2 = (eleid1 - 1)*Ky_fine*Kz_fine + (k2-1)*Ky_fine + j2;
                        [nodes_rght_face.xii2(:,:,eleid2),nodes_rght_face.eta2(:,:,eleid2),nodes_rght_face.zta2(:,:,eleid2)] = element_nodes2(nodes_rght_face.xii,nodes_rght_face.eta,nodes_rght_face.zta,i1,j1,k1,i2,j2,k2);
                    end
                end
            end
        end
    end
end

[nodes_rght_face.xxx,nodes_rght_face.yyy,nodes_rght_face.zzz] = domain_mapping(nodes_rght_face.xii2,nodes_rght_face.eta2,nodes_rght_face.zta2);

phi_rght_face = test_phi(nodes_rght_face.xxx,nodes_rght_face.yyy,nodes_rght_face.zzz);

w88_rght_face = kron(dtize_weights.wz,dtize_weights.wy);

nn_rght = +1;

Ky = Ky_coarse * Ky_fine;
Kz = Kz_coarse * Kz_fine;

for eln = 1:Ky*Kz
    boundary_dof_rght_face(:,eln) = nn_rght * reduction(phi_rght_face(:,:,eln), basis_rght_face, w88_rght_face, 1); 
end

end

