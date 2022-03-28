function [M22] = get_stiffness_matrix_3D_surfaces2(dtize_basis_surface,dtize_metric,ttl_nr_el,sz,p)

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

% M222_x1 = construct_mass_matrix2(dtize_basis_surface.yz,dtize_basis_surface.yz,dtize_metric.www,dtize_metric.g11);

% M222 = blkdiag(M222x1,M222y2,M222z3);

M22 = zeros(sz);

szt = p^2*(p+1);

M22(1:szt,1:szt,:) = M222.x1;
M22(1:szt,szt+1:2*szt,:) = M222.x2;
M22(1:szt,2*szt+1:3*szt,:) = M222.x3;

M22(szt+1:2*szt,1:szt,:) = M222.y1;
M22(szt+1:2*szt,szt+1:2*szt,:) = M222.y2;
M22(szt+1:2*szt,2*szt+1:3*szt,:) = M222.y3;

M22(2*szt+1:3*szt,1:szt,:) = M222.z1;
M22(2*szt+1:3*szt,szt+1:2*szt,:) = M222.z2;
M22(2*szt+1:3*szt,2*szt+1:3*szt,:) = M222.z3;

% M22 = [M222.x1 M222.x2 M222.x3; M222.y1 M222.y2 M222.y3; M222.z1 M222.z2 M222.z3];

end

