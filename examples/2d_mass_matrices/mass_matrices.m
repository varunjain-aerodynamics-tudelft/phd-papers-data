%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% mass_matrices.m
%
% Computes mass matrices for 2D case
%
% Written by Varun Jain - 25/03/2022
% Contact: varunjain89@gmail.com
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
close all
clear variables

%% add libraries

addpath('../../source_code')
addpath('../../source_code/MSEM')
addpath('../../source_code/invChol')
addpath('../../source_code/mass_matrix')
addpath('../../source_code/FE_functions')

%% doamin definition and mapping

c = 0.0;

xbound = [-1 1];
ybound = [-1 1];

domain.mapping  = @(xii,eta) mesh.dim_2.crazy_mesh.mapping(xbound,ybound,c,xii,eta);
domain.jacobian = @(xii,eta) mesh.dim_2.crazy_mesh.jacobian(xbound,ybound,c,xii,eta);

%% discretization

K = 1;

Kx = K;
Ky = K;

p = 3;

el_bounds.x = linspace(-1,1,Kx+1);
el_bounds.y = linspace(-1,1,Ky+1);

element.nodes = @(xii,eta,elx,ely) mesh.dim_2.element.mapping_nodes(xii,eta,elx,ely,el_bounds);

ttl_nr_el = Kx * Ky;

%% ---------- define random variables here

count.nr_volumes  = p^2;
count.nr_hor_surf = p*(p+1);
count.nr_ver_surf = p*(p+1);
count.nr_surfacs  = count.nr_hor_surf + count.nr_ver_surf;
count.nr_ttl_eln  = Kx * Ky;
count.nr_loc_lambda = 4*p;
count.nr_ttl_lambda = Kx*p*(Ky+1) + Ky*p*(Kx+1);

%% generate discretization nodes & basis

[x,w] = GLLnodes(p);
[h,e] = MimeticpolyVal(x,p,1);

dtize.basis.hx = h;
dtize.basis.hy = h;

dtize.basis.ex = e;
dtize.basis.ey = e;

dtize.weights.wx = w;
dtize.weights.wy = w;

dtize.nodes.sx = x;
dtize.nodes.sy = x;

%% generate jacobian terms on discretization nodes

[dtize.nodes.eta, dtize.nodes.xii] = meshgrid(dtize.nodes.sx,dtize.nodes.sy);

dtize.nodes.xii = dtize.nodes.xii(:)';
dtize.nodes.eta = dtize.nodes.eta(:)';

dtize.jacobian = get_jacobian_2D(p,Kx,Ky,domain.jacobian);

perm.kxx = 1;
perm.kxy = 0;
perm.kyy = 1;

dtize.metric = metric.dim_2.darcy_metric(dtize.jacobian, perm);
dtize.metric.www = kron(dtize.weights.wy,dtize.weights.wx);

%% Mass matrix M^(d-1)

dtize_basis_surface.hor = kron(dtize.basis.hy,dtize.basis.ex);
dtize_basis_surface.ver = kron(dtize.basis.ey,dtize.basis.hx);

M11 = get_stiffness_matrix_2D_surfaces(dtize_basis_surface,dtize.metric,count);

%% Mass matrix M^(d)

dtize_basis_volume = kron(dtize.basis.ey,dtize.basis.ex);

M22 = construct_mass_matrix(dtize_basis_volume,dtize_basis_volume,dtize.metric.www,dtize.metric.invg(:,:,1));