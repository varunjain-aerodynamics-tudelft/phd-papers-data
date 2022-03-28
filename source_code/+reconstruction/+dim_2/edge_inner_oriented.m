function [ux_h, uy_h ] = edge_inner_oriented(u_dof, p, rec_basis_surface, metric)

% ---------- please use this function with caution !!!!
% ---------- check your orientations and signs !!!

u1x = u_dof(1:p*(p+1));
u1y = u_dof(p*(p+1) + (1:p*(p+1)));

%% ---------- velocity x

% ---------- inner

basisX = rec_basis_surface.hor .* metric.dYdeta_g;
basisY = rec_basis_surface.ver .* metric.dYdxii_g;

ux_h = basisX .* u1x - basisY .* u1y;
ux_h = sum(ux_h);

%% ---------- velocity y

% ---------- inner

basisX = rec_basis_surface.hor .* metric.dXdeta_g;
basisY = rec_basis_surface.ver .* metric.dXdxii_g;

uy_h = -basisX .* u1x + basisY .* u1y;
uy_h = sum(uy_h);

end