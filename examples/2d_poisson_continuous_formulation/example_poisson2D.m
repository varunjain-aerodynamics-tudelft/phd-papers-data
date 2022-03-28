%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% darcy_flow_example.m
%
% Computes the pressure and velocity fields for Poisson test problem.
%
% Written by Varun Jain - 09/12/2021
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

phi_an = @(x,y) sin(2*pi*x).*sin(2*pi*y);
u_an   = @(x,y) 2*pi*cos(2*pi*x).*sin(2*pi*y);
v_an   = @(x,y) 2*pi*sin(2*pi*x).*cos(2*pi*y);
f_an   = @(x,y) -8*pi*pi*sin(2*pi*x).*sin(2*pi*y);

%% doamin definition and mapping

xbound = [-1 1];
ybound = [-1 1];
c = 0.0;

domain.mapping  = @(xii,eta) mesh.dim_2.crazy_mesh.mapping(xbound,ybound,c,xii,eta);
domain.jacobian = @(xii,eta) mesh.dim_2.crazy_mesh.jacobian(xbound,ybound,c,xii,eta,zta);

% domain_mapping = @(xi,eta) mesh.crazy_mesh.mapping(xbound, ybound, c, xi, eta);
% domain_dX_dxii = @(xi,eta) mesh.crazy_mesh.dx_dxii(xbound, c, xi, eta);
% domain_dX_deta = @(xi,eta) mesh.crazy_mesh.dx_deta(xbound, c, xi, eta);
% domain_dY_dxii = @(xi,eta) mesh.crazy_mesh.dy_dxii(ybound, c, xi, eta);
% domain_dY_deta = @(xi,eta) mesh.crazy_mesh.dy_deta(ybound, c, xi, eta);

%% discretization

K = 3;
p = 11;

element_bounds_x = linspace(-1,1,K+1);
element_bounds_y = linspace(-1,1,K+1);

element.nodes2 = @(xii,eta,zta,elx1,ely1,elz1,elx2,ely2,elz2) mesh.dim_3.element.mapping_nodes_level2(xii,eta,zta,elx1,ely1,elz1,el_bounds1,elx2,ely2,elz2,el_bounds2);

ttl_nr_el = K^2;
ttl_nr_pp = K^2*p^2;
ttl_nr_ed = 2*K*p*(K*p + 1);

local_nr_ed = 2*p*(p+1);
local_nr_pp = p^2;

%% initialize discretization nodes & basis

dtize = initialize_discretization_nodes_basis(p);

%% random variables

%% create discretization jacobian and nodes

dtize = mesh.dim_3.get_jacobian_3D_2level(dtize,Kx_coarse,Ky_coarse,Kz_coarse,Kx_fine,Ky_fine,Kz_fine,domain);

%% create discretization metric terms

%% reduction of rhs term - f

%% construct mass matrix

%% construct divergence operator / incidence matrix

%% evaluate dirichlet boundary dofs

%% construct global matrices

%% solve pressure and velocity



%% mesh mapping and derivatives

domain_mapping = @(xi,eta) mesh.crazy_mesh.mapping(xbound, ybound, c, xi, eta);
domain_dX_dxii = @(xi,eta) mesh.crazy_mesh.dx_dxii(xbound, c, xi, eta);
domain_dX_deta = @(xi,eta) mesh.crazy_mesh.dx_deta(xbound, c, xi, eta);
domain_dY_dxii = @(xi,eta) mesh.crazy_mesh.dy_dxii(ybound, c, xi, eta);
domain_dY_deta = @(xi,eta) mesh.crazy_mesh.dy_deta(ybound, c, xi, eta);

el_mapping = @(xi,eta,elx,ely) mesh.crazy_mesh.mapping_element(domain_mapping, element_bounds_x, element_bounds_y, elx, ely, xi, eta);
el_dX_dxii = @(xi,eta,elx,ely) mesh.crazy_mesh.dx_dxii_element(domain_dX_dxii, element_bounds_x, element_bounds_y, elx, ely, xi, eta);
el_dX_deta = @(xi,eta,elx,ely) mesh.crazy_mesh.dx_deta_element(domain_dX_deta, element_bounds_x, element_bounds_y, elx, ely, xi, eta);
el_dY_dxii = @(xi,eta,elx,ely) mesh.crazy_mesh.dy_dxii_element(domain_dY_dxii, element_bounds_x, element_bounds_y, elx, ely, xi, eta);
el_dY_deta = @(xi,eta,elx,ely) mesh.crazy_mesh.dy_deta_element(domain_dY_deta, element_bounds_x, element_bounds_y, elx, ely, xi, eta);

%% calculate metric terms

[xp, wp] = GLLnodes(p);
[xip, etap] = meshgrid(xp);

eval_p_dx_dxii = zeros(p+1,p+1,K^2);
eval_p_dx_deta = zeros(p+1,p+1,K^2);
eval_p_dy_dxii = zeros(p+1,p+1,K^2);
eval_p_dy_deta = zeros(p+1,p+1,K^2);

for elx = 1:K
    for ely = 1:K
        el = (elx -1)*K +ely;
        eval_p_dx_dxii(:,:,el) = el_dX_dxii(xip,etap,elx,ely);
        eval_p_dx_deta(:,:,el) = el_dX_deta(xip,etap,elx,ely);
        eval_p_dy_dxii(:,:,el) = el_dY_dxii(xip,etap,elx,ely);
        eval_p_dy_deta(:,:,el) = el_dY_deta(xip,etap,elx,ely);
    end
end

eval_p_ggg = metric.ggg(eval_p_dx_dxii, eval_p_dx_deta, eval_p_dy_dxii, eval_p_dy_deta);
eval_p_g11 = metric.g11(eval_p_dx_dxii, eval_p_dx_deta, eval_p_dy_dxii, eval_p_dy_deta, eval_p_ggg);
eval_p_g12 = metric.g12(eval_p_dx_dxii, eval_p_dx_deta, eval_p_dy_dxii, eval_p_dy_deta, eval_p_ggg);
eval_p_g22 = metric.g22(eval_p_dx_dxii, eval_p_dx_deta, eval_p_dy_dxii, eval_p_dy_deta, eval_p_ggg);

%% calculating mass matrices

M11 = zeros(local_nr_ed, local_nr_ed, ttl_nr_el);
M22 = zeros(local_nr_pp, local_nr_pp, ttl_nr_el);

for el = 1:ttl_nr_el
    M11(:,:,el) = mass_matrix.M1_3(p, eval_p_g11(:,:,el), eval_p_g12(:,:,el), eval_p_g22(:,:,el));
    M22(:,:,el) = mass_matrix.M2_2(p, eval_p_ggg(:,:,el));
end

E211 = incidence_matrix.E21(p);

M2E21 = zeros(local_nr_pp, local_nr_ed, ttl_nr_el);

for el = 1:ttl_nr_el
    M2E21(:,:,el) = M22(:,:,el) * E211;
%     M2E21(:,:,el) = E211;
end

%% calculating f-cochain

f_h = reduction.of2form_multi_element_5(f_an, p, domain_mapping, domain_dX_dxii, domain_dX_deta, domain_dY_dxii, domain_dY_deta, K, element_bounds_x, element_bounds_y);

assM2f = zeros(local_nr_pp,ttl_nr_el);

for el = 1:K^2
    assM2f(:,el) = M22(:,:,el) * f_h(:,el);
end

assM2f = assM2f(:);

%% assemble mass matrices


GM1 = gather_matrix.GM_continuous_1_form(K,p);
GM2 = gather_matrix.GM_continuous_2_form(K,p);

assM11   = AssembleMatrices2(GM1, GM1, M11);
assM2E21 = AssembleMatrices2(GM2, GM1, M2E21);

%% check stability 

try chol(assM11);
    disp('Matrix is symmetric positive definite.')
catch ME
    disp('Matrix is not symmetric positive definite')
end

d = eig(full(assM11));

isposdef = all(d > 0)

return

%% system of eq. 

LHS = [assM11 assM2E21'; assM2E21 sparse(ttl_nr_pp, ttl_nr_pp)];
RHS = [sparse(ttl_nr_ed,1); assM2f];

spy(LHS)
% export_fig('spy_cont_K3_N6_c03.pdf','-pdf','-r864','-painters','-transparent');

X_h = LHS\RHS; 

%% save X in a .mat file

%%

Q_h = X_h(1:ttl_nr_ed);
P_h = X_h(ttl_nr_ed+1:end);

q_h = zeros(local_nr_ed, ttl_nr_el);
p_h = zeros(local_nr_pp, ttl_nr_pp);

for el = 1:ttl_nr_el
    for local_edge = 1:local_nr_ed
        q_h(local_edge, el) = Q_h(GM1(local_edge, el));
    end
    for local_const = 1:local_nr_pp
        p_h(local_const, el) = P_h(GM2(local_const, el));
    end
end

%% post processing 

pf = 30;

gf = zeros(pf+1,pf+1,K^2);

[xf, wf] = GLLnodes(pf);
[xif, etaf] = meshgrid(xf);
[wfx, wfy] = meshgrid(wf);

wfxy  = wfx .* wfy;
wfxy2 = repmat(wfxy, 1, 1, K^2);

eval_pf_dx_dxii = zeros(pf+1,pf+1,K^2);
eval_pf_dx_deta = zeros(pf+1,pf+1,K^2);
eval_pf_dy_dxii = zeros(pf+1,pf+1,K^2);
eval_pf_dy_deta = zeros(pf+1,pf+1,K^2);

for elx = 1:K
    for ely = 1:K
        el = (elx -1)*K +ely;
        eval_pf_dx_dxii(:,:,el) = el_dX_dxii(xif,etaf,elx,ely);
        eval_pf_dx_deta(:,:,el) = el_dX_deta(xif,etaf,elx,ely);
        eval_pf_dy_dxii(:,:,el) = el_dY_dxii(xif,etaf,elx,ely);
        eval_pf_dy_deta(:,:,el) = el_dY_deta(xif,etaf,elx,ely);
    end
end

eval_pf_ggg = metric.ggg(eval_pf_dx_dxii, eval_pf_dx_deta, eval_pf_dy_dxii, eval_pf_dy_deta);
eval_pf_g11 = metric.g11(eval_pf_dx_dxii, eval_pf_dx_deta, eval_pf_dy_dxii, eval_pf_dy_deta, eval_pf_ggg);
eval_pf_g12 = metric.g12(eval_pf_dx_dxii, eval_pf_dx_deta, eval_pf_dy_dxii, eval_pf_dy_deta, eval_pf_ggg);
eval_pf_g22 = metric.g22(eval_pf_dx_dxii, eval_pf_dx_deta, eval_pf_dy_dxii, eval_pf_dy_deta, eval_pf_ggg);

[xf_3D, yf_3D, pf_h] = reconstruction.of2form_multi_element_v4(p_h, p, pf, el_mapping, eval_pf_ggg, K);

figure
hold on
for el = 1:K^2
    contourf(xf_3D(:,:,el), yf_3D(:,:,el), pf_h(:,:,el)')
end
colorbar

% ------- error calculation-----%

% figure
% hold on
% for el = 1:K^2
%     plot(xf_3D(:,:,el), yf_3D(:,:,el),'+')
% end



phi_an_h = phi_an(xf_3D, yf_3D);

for el = 1:K^2
    err_phi(:,:,el) = phi_an_h(:,:,el) - pf_h(:,:,el)';
end

err_phi = (err_phi).^2 .* wfxy2 .* eval_pf_ggg;
l2err_phi = sum(sum(sum(err_phi)));

l2err_phi2 = sqrt(l2err_phi)

%% reconstruction of 1-form 

half_edges = p*(p+1);

for el = 1:K^2
    qx_h(:,el) = q_h(1:half_edges,el);
    qy_h(:,el) = -q_h(half_edges+1:end,el);
end


for elx = 1:K
    for ely = 1:K
        el = (elx-1)*K + ely;

        temp = reconstruction.reconstruct1xform_2(qx_h(:,el), qy_h(:,el), p, pf, eval_pf_dx_dxii(:,:,el), eval_pf_dx_deta(:,:,el), eval_pf_dy_dxii(:,:,el), eval_pf_dy_deta(:,:,el));
        rec_qx(:,:,el) = full(temp)';
        
        temp = reconstruction.reconstruct1yform_2(qx_h(:,el), qy_h(:,el), p, pf, eval_pf_dx_dxii(:,:,el), eval_pf_dx_deta(:,:,el), eval_pf_dy_dxii(:,:,el), eval_pf_dy_deta(:,:,el));
        rec_qy(:,:,el) = full(temp)';
                
        [xf, yf] = el_mapping(xif,etaf,elx,ely);

        xf_3D2(:,:,el) = xf;
        yf_3D2(:,:,el) = yf;
    end
end

[xf, ~] = GLLnodes(pf);
[hf,ef] = MimeticpolyVal(xf,p,1);

%% for the point (-1,1) hf(:,1) and ef(:,1)

eln = 5;

%% -1,-1

dx_dxii_g_1 = eval_pf_dx_dxii(1,1,eln) / eval_pf_ggg(1,1,eln);
dx_deta_g_1 = eval_pf_dx_deta(1,1,eln) / eval_pf_ggg(1,1,eln);

dy_dxii_g_1 = eval_pf_dy_dxii(1,1,eln) / eval_pf_ggg(1,1,eln);
dy_deta_g_1 = eval_pf_dy_deta(1,1,eln) / eval_pf_ggg(1,1,eln);

for i = 1:p
    for j = 1:p+1
        edgeij = (j-1)*p + i;
        basis_hor_edge(edgeij) = ef(i,1)*hf(j,1);
    end
end

for i = 1:p+1
    for j = 1:p
        edgeij = (i-1)*p + j;
        basis_ver_edge(edgeij) = hf(i,1)*ef(j,1);
    end
end

basis_u1 = [(basis_ver_edge * dx_deta_g_1) (basis_hor_edge * dx_dxii_g_1)];
basis_v1 = [(basis_ver_edge * dy_deta_g_1) (basis_hor_edge * dy_dxii_g_1)];

u1 = basis_u1 * q_h(:,eln)
v1 = basis_v1 * q_h(:,eln)

%% exact value 

[x_1, y_1] = el_mapping(-1,-1,2,2);

u1_ex = u_an(x_1,y_1)
v1_ex = v_an(x_1,y_1)

%% 1,-1

dx_dxii_g_2 = eval_pf_dx_dxii(1,pf+1,eln) / eval_pf_ggg(1,pf+1,eln);
dx_deta_g_2 = eval_pf_dx_deta(1,pf+1,eln) / eval_pf_ggg(1,pf+1,eln);

dy_dxii_g_2 = eval_pf_dy_dxii(1,pf+1,eln) / eval_pf_ggg(1,pf+1,eln);
dy_deta_g_2 = eval_pf_dy_deta(1,pf+1,eln) / eval_pf_ggg(1,pf+1,eln);

for i = 1:p
    for j = 1:p+1
        edgeij = (j-1)*p + i;
        basis_hor_edge(edgeij) = ef(i,1)*hf(j,pf+1);
    end
end

for i = 1:p+1
    for j = 1:p
        edgeij = (i-1)*p + j;
        basis_ver_edge(edgeij) = hf(i,1)*ef(j,pf+1);
    end
end

basis_u2 = [(basis_ver_edge * dx_deta_g_2) (basis_hor_edge * dx_dxii_g_2)];
basis_v2 = [(basis_ver_edge * dy_deta_g_2) (basis_hor_edge * dy_dxii_g_2)];

u2 = basis_u2 * q_h(:,eln)
v2 = basis_v2 * q_h(:,eln)

[x_2, y_2] = el_mapping(1,-1,2,2);

u2_ex = u_an(x_2,y_2)
v2_ex = v_an(x_2,y_2)

%% 1,1

dx_dxii_g_3 = eval_pf_dx_dxii(pf+1,pf+1,eln) / eval_pf_ggg(pf+1,pf+1,eln);
dx_deta_g_3 = eval_pf_dx_deta(pf+1,pf+1,eln) / eval_pf_ggg(pf+1,pf+1,eln);

dy_dxii_g_3 = eval_pf_dy_dxii(pf+1,pf+1,eln) / eval_pf_ggg(pf+1,pf+1,eln);
dy_deta_g_3 = eval_pf_dy_deta(pf+1,pf+1,eln) / eval_pf_ggg(pf+1,pf+1,eln);

for i = 1:p
    for j = 1:p+1
        edgeij = (j-1)*p + i;
        basis_hor_edge(edgeij) = ef(i,pf+1)*hf(j,pf+1);
    end
end

for i = 1:p+1
    for j = 1:p
        edgeij = (i-1)*p + j;
        basis_ver_edge(edgeij) = hf(i,pf+1)*ef(j,pf+1);
    end
end

basis_u3 = [(basis_ver_edge * dx_deta_g_3) (basis_hor_edge * dx_dxii_g_3)];
basis_v3 = [(basis_ver_edge * dy_deta_g_3) (basis_hor_edge * dy_dxii_g_3)];

u3 = basis_u3 * q_h(:,eln)
v3 = basis_v3 * q_h(:,eln)

%% exact value 

[x_3, y_3] = el_mapping(1,1,2,2);

u3_ex = u_an(x_3,y_3)
v3_ex = v_an(x_3,y_3)


%% -1,1

dx_dxii_g_4 = eval_pf_dx_dxii(pf+1,1,eln) / eval_pf_ggg(pf+1,1,eln);
dx_deta_g_4 = eval_pf_dx_deta(pf+1,1,eln) / eval_pf_ggg(pf+1,1,eln);

dy_dxii_g_4 = eval_pf_dy_dxii(pf+1,1,eln) / eval_pf_ggg(pf+1,1,eln);
dy_deta_g_4 = eval_pf_dy_deta(pf+1,1,eln) / eval_pf_ggg(pf+1,1,eln);

for i = 1:p
    for j = 1:p+1
        edgeij = (j-1)*p + i;
        basis_hor_edge(edgeij) = ef(i,pf+1)*hf(j,1);
    end
end

for i = 1:p+1
    for j = 1:p
        edgeij = (i-1)*p + j;
        basis_ver_edge(edgeij) = hf(i,pf+1)*ef(j,1);
    end
end

basis_u4 = [(basis_ver_edge * dx_deta_g_4) (basis_hor_edge * dx_dxii_g_4)];
basis_v4 = [(basis_ver_edge * dy_deta_g_4) (basis_hor_edge * dy_dxii_g_4)];

u4 = basis_u4 * q_h(:,eln)
v4 = basis_v4 * q_h(:,eln)

[x_4, y_4] = el_mapping(-1,1,2,2);

u4_ex = u_an(x_4,y_4)
v4_ex = v_an(x_4,y_4)

%%

figure
hold on
for el = 1:K^2
    contourf(xf_3D2(:,:,el), yf_3D2(:,:,el), rec_qx(:,:,el))
end
colorbar

figure
hold on
for el = 5
    contourf(xf_3D2(:,:,el), yf_3D2(:,:,el), rec_qy(:,:,el))
end
colorbar



