%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% darcy_flow_example.m
%
% Computes the pressure and velocity fields for Darcy flow test problem in
% reference :
%   B. Wohlmuth, "A mortar finite element method using dual spaces for the
%   Lagrange multiplier", SIAM-JNA, 2000, 989-1012.
%
% Written by Varun Jain - 14/11/2020
% Contact: varunjain89@gmail.com
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% script settings

clc
close all
clear variables

%% adding library paths 

addpath('../../source_code')
addpath('../../source_code/MSEM')
addpath('../../source_code/invChol')
addpath('../../source_code/mass_matrix')
addpath('../../source_code/FE_functions')

%% test case definition 

test3.phi = @(x,y,z) x + y + z - 1.5;

test3.dpdx = @(x,y,z) 1*ones(size(x));
test3.dpdy = @(x,y,z) 1*ones(size(x));
test3.dpdz = @(x,y,z) 1*ones(size(x));

test3.d2pdx2 = @(x,y,z) zeros(size(x));
test3.d2pdy2 = @(x,y,z) zeros(size(x));
test3.d2pdz2 = @(x,y,z) zeros(size(x));

test3.d2pdxdy = @(x,y,z) zeros(size(x));
test3.d2pdxdz = @(x,y,z) zeros(size(x));
test3.d2pdydz = @(x,y,z) zeros(size(x));

test3.Kxx = @(x,y,z) x.^2 + y.^2 + 1;
test3.Kxy = @(x,y,z) zeros(size(x));
test3.Kxz = @(x,y,z) zeros(size(x));
test3.Kyy = @(x,y,z) z.^2 + 1;
test3.Kyz = @(x,y,z) sin(x .* y);
test3.Kzz = @(x,y,z) x.^2 .* y.^2 +1;

test3.dKxxdx = @(x,y,z) 2 * x;
test3.dKxydx = @(x,y,z) zeros(size(x));
test3.dKxzdx = @(x,y,z) zeros(size(x));

test3.dKxydy = @(x,y,z) zeros(size(x));
test3.dKyydy = @(x,y,z) zeros(size(x));
test3.dKyzdy = @(x,y,z) x .* cos(x .* y);

test3.dKxzdz = @(x,y,z) zeros(size(x));
test3.dKyzdz = @(x,y,z) zeros(size(x));
test3.dKzzdz = @(x,y,z) zeros(size(x));

test = test3;

test.darcy_qx = @(x,y,z) test.Kxx(x,y,z) .* test.dpdx(x,y,z) ...
                       + test.Kxy(x,y,z) .* test.dpdy(x,y,z) ...
                       + test.Kxz(x,y,z) .* test.dpdz(x,y,z);
                   
test.darcy_qy = @(x,y,z) test.Kxy(x,y,z) .* test.dpdx(x,y,z) ...
                       + test.Kyy(x,y,z) .* test.dpdy(x,y,z) ...
                       + test.Kyz(x,y,z) .* test.dpdz(x,y,z);

test.darcy_qz = @(x,y,z) test.Kxz(x,y,z) .* test.dpdx(x,y,z) ...
                       + test.Kyz(x,y,z) .* test.dpdy(x,y,z) ...
                       + test.Kzz(x,y,z) .* test.dpdz(x,y,z);
                   
test.darcy_fff = @(x,y,z) test.dKxxdx(x,y,z) .* test.dpdx(x,y,z) + test.Kxx(x,y,z) .* test.d2pdx2(x,y,z) ...
    + test.dKxydx(x,y,z) .* test.dpdy(x,y,z) + test.Kxy(x,y,z) .* test.d2pdxdy(x,y,z) ...
    + test.dKxzdx(x,y,z) .* test.dpdz(x,y,z) + test.Kxz(x,y,z) .* test.d2pdxdz(x,y,z) ...
    + test.dKxydy(x,y,z) .* test.dpdx(x,y,z) + test.Kxy(x,y,z) .* test.d2pdxdy(x,y,z) ...
    + test.dKyydy(x,y,z) .* test.dpdy(x,y,z) + test.Kyy(x,y,z) .* test.d2pdy2(x,y,z)  ...
    + test.dKyzdy(x,y,z) .* test.dpdz(x,y,z) + test.Kyz(x,y,z) .* test.d2pdydz(x,y,z) ...
    + test.dKxzdz(x,y,z) .* test.dpdx(x,y,z) + test.Kxz(x,y,z) .* test.d2pdxdz(x,y,z) ...
    + test.dKyzdz(x,y,z) .* test.dpdy(x,y,z) + test.Kyz(x,y,z) .* test.d2pdydz(x,y,z) ...
    + test.dKzzdz(x,y,z) .* test.dpdz(x,y,z) + test.Kzz(x,y,z) .* test.d2pdz2(x,y,z);

%% doamin definition and mapping

xbound = [0 1];
ybound = [0 1];
zbound = [0 1];

domain.mapping  = @(xii,eta,zta) mesh.dim_3.parallelepiped2.mapping(xbound,ybound,zbound,xii,eta,zta);
domain.jacobian = @(xii,eta,zta) mesh.dim_3.parallelepiped2.jacobian(xbound,ybound,zbound,xii,eta,zta);

%% discretization

tic

%--------- input K1 = level 1 mesh desicretization, number of elements in one direction
%--------- input K2 = level 2 mesh desicretization, number of elements in one direction

K1 = 3;
K2 = 2;

%--------- input p = order of elements

p = 3;

Kx_coarse = K1;
Ky_coarse = K1;
Kz_coarse = K1;

Kx_fine = K2;
Ky_fine = K2;
Kz_fine = K2;

el_bounds1.x = linspace(-1,1,Kx_coarse+1);
el_bounds1.y = linspace(-1,1,Ky_coarse+1);
el_bounds1.z = linspace(-1,1,Kz_coarse+1);

el_bounds2.x = linspace(-1,1,Kx_fine+1);
el_bounds2.y = linspace(-1,1,Ky_fine+1);
el_bounds2.z = linspace(-1,1,Kz_fine+1);

element.nodes2 = @(xii,eta,zta,elx1,ely1,elz1,elx2,ely2,elz2) mesh.dim_3.element.mapping_nodes_level2(xii,eta,zta,elx1,ely1,elz1,el_bounds1,elx2,ely2,elz2,el_bounds2);

%% initialize discretization nodes & basis

dtize = initialize_discretization_nodes_basis(p);

%% random variables

Kx = Kx_coarse * Kx_fine;
Ky = Ky_coarse * Ky_fine;
Kz = Kz_coarse * Kz_fine;

mesh_coarse.ttl_nr_el = Kx_coarse * Ky_coarse *Kz_coarse;

ttl_nr_el = Kx * Ky * Kz;
nr_surfc = 3*p^2*(p+1);
nr_volum = p^3;
ttl_loc_dof = nr_surfc+nr_volum;
nr_local_lambda = 6*p^2;
nr_lambda = (Kx+1)*Ky*Kz*p^2 + (Ky+1)*Kx*Kz*p^2 + (Kz+1)*Kx*Ky*p^2;

ttl_nr_mac_el = Kx_coarse*Ky_coarse*Kz_coarse;
ttl_nr_sub_el = Kx_fine*Ky_fine*Kz_fine;

ttl_nr_fine_surf = Kx_fine*p*Ky_fine*p*(Kz_fine*p+1) + Ky_fine*p*Kz_fine*p*(Kx_fine*p+1) + Kz_fine*p*Kx_fine*p*(Ky_fine*p+1);
ttl_nr_fine_voll = ttl_nr_sub_el * p^3;

q_in_coarse_el = (Kx_fine*p +1)*p^2*Ky_fine*Kz_fine... 
               + (Ky_fine*p +1)*p^2*Kx_fine*Kz_fine...
               + (Kz_fine*p +1)*p^2*Kx_fine*Ky_fine;

count_fine_nr_pressure = p^3*ttl_nr_sub_el;
nr_boundary_lambda = 2*p^2*(Ky_fine*Kz_fine + Kx_fine*Kz_fine + Kx_fine*Ky_fine);

nr_int_lambda = Kx_coarse*Kx_fine*p*Ky_coarse*Ky_fine*p*(Kz_coarse-1) + Ky_fine*Ky_coarse*p*Kz_fine*Kz_coarse*p*(Kx_coarse-1) + Kx_fine*Kx_coarse*p*Kz_fine*Kz_coarse*p*(Ky_coarse-1);
nr_out_lambda = 2*Kx_fine*p*Ky_fine*p * Kx_coarse*Ky_coarse + 2*Ky_fine*p*Kz_fine*p * Ky_coarse*Kz_coarse + 2*Kx_fine*p*Kz_fine*p  * Kx_coarse*Kz_coarse;

ttl_nr_sub_surfc = p^2*Kx_fine*Ky_fine*(Kz_fine*p+1) + p^2*Ky_fine*Kz_fine*(Kx_fine*p+1) + p^2*Kx_fine*Kz_fine*(Ky_fine*p+1);

%% create discretization jacobian and nodes

dtize = mesh.dim_3.get_jacobian_3D_2level(dtize,Kx_coarse,Ky_coarse,Kz_coarse,Kx_fine,Ky_fine,Kz_fine,domain);

%% create discretization metric terms

%--------- evaluate permiability tensor at discretization nodes

perm.kxx = test.Kxx(dtize.nodes.x, dtize.nodes.y, dtize.nodes.z);
perm.kxy = test.Kxy(dtize.nodes.x, dtize.nodes.y, dtize.nodes.z);
perm.kxz = test.Kxz(dtize.nodes.x, dtize.nodes.y, dtize.nodes.z);
perm.kyy = test.Kyy(dtize.nodes.x, dtize.nodes.y, dtize.nodes.z);
perm.kyz = test.Kyz(dtize.nodes.x, dtize.nodes.y, dtize.nodes.z);
perm.kzz = test.Kzz(dtize.nodes.x, dtize.nodes.y, dtize.nodes.z);

%--------- evaluate metric tensor at discretization nodes

dtize_metric = metric.dim_3.get_metric_3D_darcy(dtize.jacobian,perm);
dtize_metric.www = kron(kron(dtize.weights.wz,dtize.weights.wy),dtize.weights.wx);

dtize_metric_red.ggg = metric.dim_3.determinant_jacobian(dtize.jacobian);
dtize_metric_red.invg = 1 ./dtize_metric_red.ggg;
dtize_metric_red = metric.dim_3.poisson_metric(dtize.jacobian, dtize_metric_red.invg);
dtize_metric_red.ggg = metric.dim_3.determinant_jacobian(dtize.jacobian);
dtize_metric_red.invg = 1 ./dtize_metric_red.ggg;
dtize_metric_red.www = kron(kron(dtize.weights.wz,dtize.weights.wy),dtize.weights.wx);

clear perm

%% reduction of rhs term - f

basis_volume = kron(kron(dtize.basis.ez,dtize.basis.ey),dtize.basis.ex);

M333 = zeros(nr_volum, nr_volum, ttl_nr_el);
for eln = 1:ttl_nr_el
    M333(:,:,eln) = construct_mass_matrix(basis_volume,basis_volume,dtize.www,dtize_metric.invg(:,:,eln));
end

ff_int = test.darcy_fff(dtize.nodes.x,dtize.nodes.y,dtize.nodes.z);
dof_f = zeros(nr_volum,1,ttl_nr_el);
for eln = 1:ttl_nr_el
    dof_f(:,1,eln) = chol_sol(M333(:,:,eln),reduction(ff_int(:,:,eln), basis_volume, dtize.www,1));
end

% clear M333

%% construct mass matrix

dtize.basis.surface.yz = kron(kron(dtize.basis.ez,dtize.basis.ey),dtize.basis.hx);
dtize.basis.surface.zx = kron(kron(dtize.basis.ez,dtize.basis.hy),dtize.basis.ex);
dtize.basis.surface.xy = kron(kron(dtize.basis.hz,dtize.basis.ey),dtize.basis.ex);

sz = [nr_surfc nr_surfc ttl_nr_el];
M22 = get_stiffness_matrix_3D_surfaces2(dtize.basis.surface,dtize_metric,ttl_nr_el,sz,p);
M22_red = get_stiffness_matrix_3D_surfaces2(dtize.basis.surface,dtize_metric_red,ttl_nr_el,sz,p);

clear dtize_metric

%% construct divergence operator / incidence matrix

E32 = incidence_matrix.dim_3.discrete_divergence(p);

%% evaluate dirichlet boundary dofs

local_bndry_face = mesh.dim_3.element_boundary_faces_map(p);
[boundary_global_ID] = mesh.dim_3.get_hybrid_1_boundary_element_ID(Kx_coarse,Ky_coarse,Kz_coarse,Kx_fine,Ky_fine,Kz_fine);

dof_dirichlet_bc = zeros(nr_surfc,1,ttl_nr_el);
dof_dirichlet_bc(local_bndry_face.left,1,boundary_global_ID.left) = get_boundary_dof.hybrid_1.left(dtize.basis, dtize.nodes, Kx_coarse, Ky_coarse, Kz_coarse, Kx_fine, Ky_fine, Kz_fine, test.phi, dtize.weights, element.nodes2, domain.mapping);
dof_dirichlet_bc(local_bndry_face.rght,1,boundary_global_ID.rght) = get_boundary_dof.hybrid_1.rght(dtize.basis, dtize.nodes, Kx_coarse, Ky_coarse, Kz_coarse, Kx_fine, Ky_fine, Kz_fine, test.phi, dtize.weights, element.nodes2, domain.mapping);

%% construct global matrices

% --------- generate gather matrices ---------%

gather_qqq = gather_matrix.dim_3.continuous_surfaces(Kx_fine, Ky_fine, Kz_fine, p);

gather_ppp = zeros(nr_volum,ttl_nr_sub_el);
for eln = 1:ttl_nr_sub_el
    gather_ppp(1:nr_volum,eln) = (eln-1)*nr_volum + (1:nr_volum);
end

%--------- assemble divergence operator ---------%

E33 = repmat(full(E32), [1 1 ttl_nr_sub_el]);
coll_E32 = AssembleMatrices.AssembleMatrices2(gather_ppp, gather_qqq, E33);

%--------- assembled rhs term - f ---------%

F_h = zeros(nr_volum*ttl_nr_sub_el,ttl_nr_mac_el);

for eln = 1:mesh_coarse.ttl_nr_el
    temp = dof_f(:,1,(eln-1)*ttl_nr_sub_el+1:eln*ttl_nr_sub_el);
    F_h(:,eln) = temp(:);
end

%--------- assemble dirichlet boundary conditions ---------%

B_D = zeros(q_in_coarse_el,ttl_nr_mac_el);

for eln = 1:mesh_coarse.ttl_nr_el
    temp = dof_dirichlet_bc(:,:,(eln-1)*ttl_nr_sub_el + 1:eln*ttl_nr_sub_el);
    B_D(:,eln) = AssembleMatrices.AssembleMatrices2(gather_qqq,ones(1,ttl_nr_sub_el),temp);
end

clear dof_dirichlet_bc

%--------- construct local Neumann boundary

basis_surfac_yz = kron(kron(dtize.basis.ez,dtize.basis.ey),dtize.basis.hx);
basis_surfac_zx = kron(kron(dtize.basis.ez,dtize.basis.hy),dtize.basis.ex);
basis_surfac_xy = kron(kron(dtize.basis.hz,dtize.basis.ey),dtize.basis.ex);

qqx_int = test.darcy_qx(dtize.nodes.x,dtize.nodes.y,dtize.nodes.z);
qqy_int = test.darcy_qy(dtize.nodes.x,dtize.nodes.y,dtize.nodes.z);
qqz_int = test.darcy_qz(dtize.nodes.x,dtize.nodes.y,dtize.nodes.z);
    
red.metric.x = qqx_int .* dtize.jacobian.dXdxii + qqy_int .* dtize.jacobian.dYdxii + qqz_int .* dtize.jacobian.dZdxii;
red.metric.y = qqx_int .* dtize.jacobian.dXdeta + qqy_int .* dtize.jacobian.dYdeta + qqz_int .* dtize.jacobian.dZdeta;
red.metric.z = qqx_int .* dtize.jacobian.dXdzta + qqy_int .* dtize.jacobian.dYdzta + qqz_int .* dtize.jacobian.dZdzta;

dof_uu = zeros(nr_surfc,ttl_nr_fine_voll);
for eln = 1:Kx*Ky*Kz
    dof_xx = reduction(red.metric.x(:,:,eln), basis_surfac_yz, dtize.www, 1);
    dof_yy = reduction(red.metric.y(:,:,eln), basis_surfac_zx, dtize.www, 1);
    dof_zz = reduction(red.metric.z(:,:,eln), basis_surfac_xy, dtize.www, 1);
    
    dof_temp = [dof_xx; dof_yy; dof_zz];
    dof_uu(:,eln) = M22_red(:,:,eln) \ dof_temp;
end

%--------- gather dof from local to global

for eln1 = 1:mesh_coarse.ttl_nr_el
   eleID = (eln1 - 1)*Kx_fine*Ky_fine*Kz_fine + (1:Kx_fine*Ky_fine*Kz_fine);
   temp_dof_topp = reshape(dof_uu(:,eleID), [nr_surfc 1 Kx_fine*Ky_fine*Kz_fine]);
   dof_X.U(:,eln1) = AssembleMatrices.AssembleMatrices2_full_no_add(gather_qqq,ones(1,Kx_fine*Ky_fine*Kz_fine),temp_dof_topp);
end

bounndary_eleID = mesh.dim_3.get_boundary_element_ID(Kx_coarse,Ky_coarse, Kz_coarse);

local_N = incidence_matrix.dim_3.discrete_unit_normal_multiple_element_cont_new(p,Kx_fine,Ky_fine,Kz_fine,local_bndry_face.rght,local_bndry_face.left,local_bndry_face.topp,local_bndry_face.bttm,local_bndry_face.frnt,local_bndry_face.back);
local_N_topp = incidence_matrix.dim_3.unit_normal_top_boundary(p,Kx_fine,Ky_fine,Kz_fine,local_bndry_face.rght,local_bndry_face.left,local_bndry_face.topp,local_bndry_face.bttm,local_bndry_face.frnt,local_bndry_face.back);
local_N_bttm = incidence_matrix.dim_3.unit_normal_bot_boundary(p,Kx_fine,Ky_fine,Kz_fine,local_bndry_face.rght,local_bndry_face.left,local_bndry_face.topp,local_bndry_face.bttm,local_bndry_face.frnt,local_bndry_face.back);
local_N_frnt = incidence_matrix.dim_3.unit_normal_frnt_boundary(p,Kx_fine,Ky_fine,Kz_fine,local_bndry_face.rght,local_bndry_face.left,local_bndry_face.topp,local_bndry_face.bttm,local_bndry_face.frnt,local_bndry_face.back);
local_N_back = incidence_matrix.dim_3.unit_normal_back_boundary(p,Kx_fine,Ky_fine,Kz_fine,local_bndry_face.rght,local_bndry_face.left,local_bndry_face.topp,local_bndry_face.bttm,local_bndry_face.frnt,local_bndry_face.back);

B_N_topp = zeros(nr_boundary_lambda,mesh_coarse.ttl_nr_el);
B_N_bttm = zeros(nr_boundary_lambda,mesh_coarse.ttl_nr_el);
B_N_frnt = zeros(nr_boundary_lambda,mesh_coarse.ttl_nr_el);
B_N_back = zeros(nr_boundary_lambda,mesh_coarse.ttl_nr_el);

B_N_topp(:,bounndary_eleID.topp) = local_N * (local_N_topp' * local_N_topp) *dof_X.U(:,bounndary_eleID.topp);
B_N_bttm(:,bounndary_eleID.bttm) = local_N * (local_N_bttm' * local_N_bttm) *dof_X.U(:,bounndary_eleID.bttm);
B_N_frnt(:,bounndary_eleID.frnt) = local_N * (local_N_frnt' * local_N_frnt) *dof_X.U(:,bounndary_eleID.frnt);
B_N_back(:,bounndary_eleID.back) = local_N * (local_N_back' * local_N_back) *dof_X.U(:,bounndary_eleID.back);

B_N = B_N_topp + B_N_bttm + B_N_frnt + B_N_back;

%% --------- evaluate inverse of M22 ---------%

inv_M22_coarse = zeros(q_in_coarse_el,q_in_coarse_el,mesh_coarse.ttl_nr_el);
inv_EME = zeros(count_fine_nr_pressure,count_fine_nr_pressure,mesh_coarse.ttl_nr_el);

for eln = 1:mesh_coarse.ttl_nr_el
    temp = AssembleMatrices.AssembleMatrices2_full(gather_qqq, gather_qqq, M22(:,:,(eln-1)*ttl_nr_sub_el + 1: eln * ttl_nr_sub_el));
    inv_M22_coarse(:,:,eln) = invChol_mex(temp);
    inv_EME(:,:,eln) = invChol_mex(coll_E32 * inv_M22_coarse(:,:,eln) * coll_E32');
end

clear M22_coarse

%% --------- construct lambda system ---------%

lam_lhs = zeros(nr_boundary_lambda,nr_boundary_lambda,mesh_coarse.ttl_nr_el);
lam_rhs = zeros(nr_boundary_lambda,1,mesh_coarse.ttl_nr_el);

for eln = 1:mesh_coarse.ttl_nr_el
    %--------- construct LHS
    lam_1 = local_N * inv_M22_coarse(:,:,eln) * local_N';
    lam_2 = coll_E32 * inv_M22_coarse(:,:,eln) * local_N';

    lam_lhs(:,:,eln) = lam_1 - lam_2' * inv_EME(:,:,eln) * lam_2;
    
    %--------- construct RHS
    rhs_lam1 = coll_E32' * (inv_EME(:,:,eln) * F_h(:,eln));
    rhs_lam2 = -coll_E32' * (inv_EME(:,:,eln) * (coll_E32 * (inv_M22_coarse(:,:,eln) * B_D(:,eln))));

    lam_rhs(:,1,eln) = -B_N(:,eln) + local_N * (inv_M22_coarse(:,:,eln) * (rhs_lam1 + rhs_lam2 + B_D(:,eln)));
end

clear lam_1 lam_2 rhs_lam1 rhs_lam2  

%--------- assemble lambda system

gather_lambda = gather_matrix.dim_3.lambda_pressure_multiple_element3(Kx_coarse,Ky_coarse,Kz_coarse,Kx_fine,Ky_fine,Kz_fine,p);

LHS_lambda = AssembleMatrices.AssembleMatrices2(gather_lambda,gather_lambda,lam_lhs);
RHS_lambda = full(AssembleMatrices.AssembleMatrices2(gather_lambda, ones(1,ttl_nr_mac_el), lam_rhs));

clear lam_lhs lam_rhs

% --------- remove Dirichlet Bc's

rght_face_lambda_ID = nr_int_lambda + (1:Ky_fine*p*Kz_fine*p * Ky_coarse*Kz_coarse);
left_face_lambda_ID = nr_int_lambda + Ky_fine*p*Kz_fine*p * Ky_coarse*Kz_coarse + (1:Ky_fine*p*Kz_fine*p * Ky_coarse*Kz_coarse);

dirichlet_boundary_ID = [rght_face_lambda_ID left_face_lambda_ID];

LHS_lambda(dirichlet_boundary_ID,:) = [];
LHS_lambda(:,dirichlet_boundary_ID) = [];
RHS_lambda(dirichlet_boundary_ID,:) = [];

%% solve lambda 

%--------- solve lambda system using cholesky decomposition

dLHS_lambda = decomposition(LHS_lambda,'chol','lower');
X.L = dLHS_lambda\RHS_lambda;

% --------- add Dirichlet boundaries back 

nr_rght_face = Ky_fine*p*Kz_fine*p * Ky_coarse*Kz_coarse;
nr_topp_face = Kx_fine*p*Kz_fine*p * Kx_coarse*Kz_coarse;
nr_frnt_face = Kx_fine*p*Ky_fine*p * Kx_coarse*Ky_coarse;

X.L = [X.L(1:nr_int_lambda); zeros(2*nr_rght_face,1); X.L(nr_int_lambda + (1:2*nr_topp_face+2*nr_frnt_face))];
X.L = X.L(gather_lambda);

%% solve for pressure and velocity

for eln = 1:mesh_coarse.ttl_nr_el
    X_p = coll_E32 * (inv_M22_coarse(:,:,eln) * (B_D(:,eln) - local_N' *X.L(:,eln)));
    X.P(:,eln) = inv_EME(:,:,eln) * (-F_h(:,eln) + X_p);
    
    X.U(:,eln) = inv_M22_coarse(:,:,eln) * (B_D(:,eln) - coll_E32'*X.P(:,eln) - local_N' * X.L(:,eln));
end

time = toc;

check2 = X.U - dof_X.U;

%% post processing

pf = p+5;

%% create reconstruction nodes & basis

[xf,wf] = GLLnodes(pf);
[hf,ef] = MimeticpolyVal(xf,p,1);

rec.basis.hfx = hf;
rec.basis.hfy = hf;
rec.basis.hfz = hf;

rec.basis.efx = ef;
rec.basis.efy = ef;
rec.basis.efz = ef;

rec.weights.wfx = wf;
rec.weights.wfy = wf;
rec.weights.wfz = wf;

rec.nodes.sfx = xf;
rec.nodes.sfy = xf;
rec.nodes.sfz = xf;

rec.basis.hx = hf;
rec.basis.hy = hf;
rec.basis.hz = hf;

rec.basis.ex = ef;
rec.basis.ey = ef;
rec.basis.ez = ef;

rec.weights.wx = wf;
rec.weights.wy = wf;
rec.weights.wz = wf;

rec.nodes.sx = xf;
rec.nodes.sy = xf;
rec.nodes.sz = xf;

rec_basis_volume = kron(kron(rec.basis.efz,rec.basis.efy),rec.basis.efx);

%------- create reconstruction jacobian and nodes

rec = mesh.dim_3.get_jacobian_3D_2level(rec,Kx_coarse,Ky_coarse,Kz_coarse,Kx_fine,Ky_fine,Kz_fine,domain);
rec.metric.ggg = metric.dim_3.determinant_jacobian(rec.jacobian);
rec.wwwf = kron(kron(rec.weights.wfz,rec.weights.wfy),rec.weights.wfx);

%----------- interp errror in F :: f - f^h

for eln = 1:ttl_nr_el
    rec.fff(:,:,eln) = reconstruction.dim_3.volume(dof_f(:,1,eln),rec_basis_volume,pf,rec.metric.ggg(:,:,eln))';
end

ext_fff = test.darcy_fff(rec.nodes.x,rec.nodes.y,rec.nodes.z);
error.fff = error_processor_v3(ext_fff,rec.fff,rec.wwwf,rec.metric.ggg);

disp('error in l2 norm - projection of f')
error.fff.sqrt

%---------- interp errror in f - div_q^h

xx.uu = zeros(nr_surfc,ttl_nr_el);
xx.pp = zeros(nr_volum,ttl_nr_el);

for eln = 1:ttl_nr_mac_el
    temp1 = X.U(:,eln);
    xx.uu(:,(eln-1)*ttl_nr_sub_el+1:eln*ttl_nr_sub_el) = temp1(gather_qqq);
    
    temp2 = X.P(:,eln);
    xx.pp(:,(eln-1)*ttl_nr_sub_el+1:eln*ttl_nr_sub_el) = temp2(gather_ppp);
end

for eln = 1:ttl_nr_mac_el
    temp1 = dof_X.U(:,eln);
    xx2.uu(:,(eln-1)*ttl_nr_sub_el+1:eln*ttl_nr_sub_el) = temp1(gather_qqq);
end

% !!!!!!!!!!! check 
check1 = xx2.uu - dof_uu;

div_q = zeros(nr_volum,ttl_nr_el);

for eln = 1:ttl_nr_el
    div_q(:,eln) = E32 * xx.uu(:,eln);
end

for eln = 1:ttl_nr_el
    rec.div_q(:,:,eln) = reconstruction.dim_3.volume(div_q(:,eln),rec_basis_volume,pf,rec.metric.ggg(:,:,eln))';
end

error.div_q = error_processor_v3(ext_fff,rec.div_q,rec.wwwf,rec.metric.ggg);

disp('error in l2 norm - divergence of velocity')
error.div_q.sqrt

%----------- constraint error ::: f^h - div q^h

rec.f_h_div_q_h = rec.fff - rec.div_q;
error.f_h_div_q_h = error_processor_v3(0,rec.f_h_div_q_h,rec.wwwf,rec.metric.ggg);
error.f_h_div_q_h.sqrt;

%------------- l2 error in p

%----- change dual dof to primal 

for eln = 1:ttl_nr_el
    xx.pp(:,eln) = chol_sol(M333(:,:,eln),xx.pp(:,eln));
    rec.ppp(:,:,eln) = reconstruction.dim_3.volume(xx.pp(:,eln),rec_basis_volume,pf,rec.metric.ggg(:,:,eln))';
end

ext_ppp = test.phi(rec.nodes.x,rec.nodes.y,rec.nodes.z);
error.ppp = error_processor_v3(ext_ppp,rec.ppp,rec.wwwf,rec.metric.ggg);
disp('error in l2 norm - pressure')
error.ppp.sqrt

X.L2 = X.L - local_N * B_D;

M22_red_coarse = zeros(ttl_nr_fine_surf,ttl_nr_fine_surf,mesh_coarse.ttl_nr_el);
for eln = 1:mesh_coarse.ttl_nr_el
    M22_red_coarse(:,:,eln) = AssembleMatrices.AssembleMatrices2_full(gather_qqq, gather_qqq, M22_red(:,:,(eln-1)*ttl_nr_sub_el + 1: eln * ttl_nr_sub_el));
end

grad_p = zeros(ttl_nr_fine_surf,mesh_coarse.ttl_nr_el);
for eln1 = 1:mesh_coarse.ttl_nr_el
    grad_p(:,eln1) = M22_red_coarse(:,:,eln1)\(- coll_E32' * X.P(:,eln1) - local_N' * X.L2(:,eln1));
end

for eln = 1:ttl_nr_mac_el
    temp1 = grad_p(:,eln);
    xx_grad_p.uu(:,(eln-1)*ttl_nr_sub_el+1:eln*ttl_nr_sub_el) = temp1(gather_qqq);
end

%------------ H div error in u

rec.basis.surfac.xy = kron(kron(rec.basis.hfz,rec.basis.efy),rec.basis.efx);
rec.basis.surfac.yz = kron(kron(rec.basis.efz,rec.basis.efy),rec.basis.hfx);
rec.basis.surfac.zx = kron(kron(rec.basis.efz,rec.basis.hfy),rec.basis.efx);

rec.metric.invg = 1./rec.metric.ggg;

for eln = 1:ttl_nr_el
    
    local_jac.dXdxii = rec.jacobian.dXdxii(:,:,eln);
    local_jac.dXdeta = rec.jacobian.dXdeta(:,:,eln);
    local_jac.dXdzta = rec.jacobian.dXdzta(:,:,eln);
    local_jac.dYdxii = rec.jacobian.dYdxii(:,:,eln);
    local_jac.dYdeta = rec.jacobian.dYdeta(:,:,eln);
    local_jac.dYdzta = rec.jacobian.dYdzta(:,:,eln);
    local_jac.dZdxii = rec.jacobian.dZdxii(:,:,eln);
    local_jac.dZdeta = rec.jacobian.dZdeta(:,:,eln);
    local_jac.dZdzta = rec.jacobian.dZdzta(:,:,eln);
    
    [rec.uu(:,:,eln), rec.vv(:,:,eln), rec.ww(:,:,eln)] = reconstruction.dim_3.surface2(xx.uu(:,eln), local_jac, p,rec.metric.invg(:,:,eln),rec.basis.surfac);
    [rec_grad_p.uu(:,:,eln), rec_grad_p.vv(:,:,eln), rec_grad_p.ww(:,:,eln)] = reconstruction.dim_3.surface2(xx_grad_p.uu(:,eln), local_jac, p,rec.metric.invg(:,:,eln),rec.basis.surfac);
end

ext_uu = test.darcy_qx(rec.nodes.x,rec.nodes.y,rec.nodes.z);
ext_vv = test.darcy_qy(rec.nodes.x,rec.nodes.y,rec.nodes.z);
ext_ww = test.darcy_qz(rec.nodes.x,rec.nodes.y,rec.nodes.z);

ext_grad_p_uu = test.dpdx(rec.nodes.x,rec.nodes.y,rec.nodes.z);
ext_grad_p_vv = test.dpdy(rec.nodes.x,rec.nodes.y,rec.nodes.z);
ext_grad_p_ww = test.dpdz(rec.nodes.x,rec.nodes.y,rec.nodes.z);

error.uuu = error_processor_v3(ext_uu, rec.uu, rec.wwwf, rec.metric.ggg);
error.vvv = error_processor_v3(ext_vv, rec.vv, rec.wwwf, rec.metric.ggg);
error.www = error_processor_v3(ext_ww, rec.ww, rec.wwwf, rec.metric.ggg);

error_grad_p.uuu = error_processor_v3(ext_grad_p_uu, rec_grad_p.uu, rec.wwwf, rec.metric.ggg);
error_grad_p.vvv = error_processor_v3(ext_grad_p_vv, rec_grad_p.vv, rec.wwwf, rec.metric.ggg);
error_grad_p.www = error_processor_v3(ext_grad_p_ww, rec_grad_p.ww, rec.wwwf, rec.metric.ggg);

error.uuu.sqrt;
error.vvv.sqrt;
error.www.sqrt;

error_grad_p.uuu.sqrt;
error_grad_p.vvv.sqrt;
error_grad_p.www.sqrt;

error.q_Hdiv.sqrt = sqrt(error.uuu.sqre + error.vvv.sqre + error.www.sqre + error.div_q.sqre);

disp('error in Hdiv norm - velocity')
error.q_Hdiv.sqrt

%% plots

[boundary_global_ID] = mesh.dim_3.get_hybrid_1_boundary_element_ID(Kx_coarse,Ky_coarse,Kz_coarse,Kx_fine,Ky_fine,Kz_fine);

mesh.dim_3.plot_3D_hybrid_0(rec.ppp,rec.nodes,Kx,Ky,Kz,boundary_global_ID,pf)
title('pressure field')
mesh.dim_3.plot_3D_hybrid_0(rec.uu,rec.nodes,Kx,Ky,Kz,boundary_global_ID,pf)
title('velocity field in x-direction')
mesh.dim_3.plot_3D_hybrid_0(rec.vv,rec.nodes,Kx,Ky,Kz,boundary_global_ID,pf)
title('velocity field in y-direction')
mesh.dim_3.plot_3D_hybrid_0(rec.ww,rec.nodes,Kx,Ky,Kz,boundary_global_ID,pf)
title('velocity field in z-direction')