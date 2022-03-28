function [M11] = get_stiffness_matrix_2D_surfaces(dtize_basis_surface,dtize_metric,count)

% count_hor_surf  = count.nr_hor_surf;
% count_ver_surf  = count.nr_ver_surf;
count_ttl_nr_el = count.nr_ttl_eln;

for eln = 1:count_ttl_nr_el
    M111.x1(:,:,eln) = construct_mass_matrix(dtize_basis_surface.hor,dtize_basis_surface.hor,dtize_metric.www,dtize_metric.g11(:,:,eln));
    M111.x2(:,:,eln) = construct_mass_matrix(dtize_basis_surface.hor,dtize_basis_surface.ver,dtize_metric.www,dtize_metric.g12(:,:,eln));

    M111.y1(:,:,eln) = construct_mass_matrix(dtize_basis_surface.ver,dtize_basis_surface.hor,dtize_metric.www,dtize_metric.g12(:,:,eln));
    M111.y2(:,:,eln) = construct_mass_matrix(dtize_basis_surface.ver,dtize_basis_surface.ver,dtize_metric.www,dtize_metric.g22(:,:,eln));    
end

M11 = [M111.x1 M111.x2; M111.y1 M111.y2];

% count_nr_surf = count_hor_surf + count_ver_surf;
% 
% M11 = zeros(count_nr_surf,count_nr_surf,count_ttl_nr_el);
% 
% M11(1:count_hor_surf,1:count_hor_surf,:) = M111.x1;
% M11(1:count_hor_surf,count_hor_surf+1:count_nr_surf,:) = M111.x2;
% 
% M11(count_hor_surf+1:count_nr_surf,1:count_hor_surf,:) = M111.y1;
% M11(count_hor_surf+1:count_nr_surf,1:count_hor_surf,:) = M111.y2;

end

