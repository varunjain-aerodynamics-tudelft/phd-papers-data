function [M22] = get_stiffness_matrix_3D_surfaces(dtize_basis_surface,dtize_metric,ttl_nr_el)

for eln = 1:ttl_nr_el
    M222.x1(:,:,eln) = construct_mass_matrix(dtize_basis_surface.yz,dtize_basis_surface.yz,dtize_metric.www,dtize_metric.g11(:,:,eln));
    M222.x2(:,:,eln) = construct_mass_matrix(dtize_basis_surface.yz,dtize_basis_surface.zx,dtize_metric.www,dtize_metric.g12(:,:,eln));
    M222.x3(:,:,eln) = construct_mass_matrix(dtize_basis_surface.yz,dtize_basis_surface.xy,dtize_metric.www,dtize_metric.g13(:,:,eln));
    
    M222.y1(:,:,eln) = construct_mass_matrix(dtize_basis_surface.zx,dtize_basis_surface.yz,dtize_metric.www,dtize_metric.g12(:,:,eln));
    M222.y2(:,:,eln) = construct_mass_matrix(dtize_basis_surface.zx,dtize_basis_surface.zx,dtize_metric.www,dtize_metric.g22(:,:,eln));
    M222.y3(:,:,eln) = construct_mass_matrix(dtize_basis_surface.zx,dtize_basis_surface.xy,dtize_metric.www,dtize_metric.g23(:,:,eln));
    
    M222.z1(:,:,eln) = construct_mass_matrix(dtize_basis_surface.xy,dtize_basis_surface.yz,dtize_metric.www,dtize_metric.g13(:,:,eln));
    M222.z2(:,:,eln) = construct_mass_matrix(dtize_basis_surface.xy,dtize_basis_surface.zx,dtize_metric.www,dtize_metric.g23(:,:,eln));
    M222.z3(:,:,eln) = construct_mass_matrix(dtize_basis_surface.xy,dtize_basis_surface.xy,dtize_metric.www,dtize_metric.g33(:,:,eln));
end

% M222 = blkdiag(M222x1,M222y2,M222z3);

M22 = [M222.x1 M222.x2 M222.x3; M222.y1 M222.y2 M222.y3; M222.z1 M222.z2 M222.z3];

end

