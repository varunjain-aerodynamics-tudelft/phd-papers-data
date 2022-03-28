function [ux_h, uy_h ] = edge_outer_oriented(u_dof, p, rec_basis_surface, metric)

u1x = -u_dof(1:p*(p+1));
u1y = u_dof(p*(p+1) + (1:p*(p+1)));

%% ---------- velocity x

% ---------- inner

basisX = rec_basis_surface.hor .* metric.dYdeta_g;
basisY = rec_basis_surface.ver .* metric.dYdxii_g;

uy_h = basisX .* u1x - basisY .* u1y;
uy_h = -sum(uy_h);

% ---------- outer

% basisX = rec_basis_surface.hor .* metric.dXdxii_g;
% basisY = rec_basis_surface.ver .* metric.dXdeta_g;
% 
% ux_h = - basisX .* u1x + basisY .* u1y;
% ux_h = sum(ux_h);

%% ---------- velocity y

% ---------- inner

basisX = rec_basis_surface.hor .* metric.dXdeta_g;
basisY = rec_basis_surface.ver .* metric.dXdxii_g;

ux_h = -basisX .* u1x + basisY .* u1y;
ux_h = sum(ux_h);

% ---------- outer

% basisX = rec_basis_surface.hor .* metric.dYdxii_g;
% basisY = rec_basis_surface.ver .* metric.dYdeta_g;
% 
% uy_h = - basisX .* u1x - basisY .* u1y;
% uy_h = -sum(uy_h);

end