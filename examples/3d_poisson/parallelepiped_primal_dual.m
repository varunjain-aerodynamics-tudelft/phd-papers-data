%% script settings

disp('clearing variables')
disp('...')

clc
close all
clear variables

%% adding library paths 

disp('adding library paths')
disp('...')

addpath('../../source_code/')
addpath('../../source_code/MSEM')
addpath('../../source_code/mass_matrix')
addpath('../../source_code/FE_functions')
% addpath('../../../../export_fig')

%% test case

disp('define test case')
disp('...')

% test 3

test3.phi = @(x,y,z) sin(2*pi*x) .* sin(2*pi*y) .* sin(2*pi*z);

test3.dpdx = @(x,y,z) 2*pi * cos(2*pi*x) .* sin(2*pi*y) .* sin(2*pi*z);
test3.dpdy = @(x,y,z) 2*pi * sin(2*pi*x) .* cos(2*pi*y) .* sin(2*pi*z);
test3.dpdz = @(x,y,z) 2*pi * sin(2*pi*x) .* sin(2*pi*y) .* cos(2*pi*z);

test3.d2pdx2 = @(x,y,z) -4*pi^2 * sin(2*pi*x) .* sin(2*pi*y) .* sin(2*pi*z);
test3.d2pdy2 = @(x,y,z) -4*pi^2 * sin(2*pi*x) .* sin(2*pi*y) .* sin(2*pi*z);
test3.d2pdz2 = @(x,y,z) -4*pi^2 * sin(2*pi*x) .* sin(2*pi*y) .* sin(2*pi*z);

test3.d2pdxdy = @(x,y,z) zeros(size(x));
test3.d2pdxdz = @(x,y,z) zeros(size(x));
test3.d2pdydz = @(x,y,z) zeros(size(x));

%% 

test3.Kxx = @(x,y,z) 1;
test3.Kxy = @(x,y,z) 0;
test3.Kxz = @(x,y,z) 0;
test3.Kyy = @(x,y,z) 1;
test3.Kyz = @(x,y,z) 0;
test3.Kzz = @(x,y,z) 1;

test3.dKxxdx = @(x,y,z) 0;
test3.dKxydx = @(x,y,z) 0;
test3.dKxzdx = @(x,y,z) 0;

test3.dKxydy = @(x,y,z) 0;
test3.dKyydy = @(x,y,z) 0;
test3.dKyzdy = @(x,y,z) 0;

test3.dKxzdz = @(x,y,z) 0;
test3.dKyzdz = @(x,y,z) 0;
test3.dKzzdz = @(x,y,z) 0;

%% 

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

disp('define domain mapping')
disp('...')

xbound = [0 1];
ybound = [0 1];
zbound = [0 1];

% c = 0.1;
% 
% domain.mapping  = @(xii,eta,zta) mesh.dim_3.crazy_cube.mapping(xbound,ybound,zbound,c,xii,eta,zta);
% domain.jacobian = @(xii,eta,zta) mesh.dim_3.crazy_cube.jacobian(xbound,ybound,zbound,c,xii,eta,zta);

domain.mapping  = @(xii,eta,zta) mesh.dim_3.parallelepiped2.mapping(xbound,ybound,zbound,xii,eta,zta);
domain.jacobian = @(xii,eta,zta) mesh.dim_3.parallelepiped2.jacobian(xbound,ybound,zbound,xii,eta,zta);

%% discretization
tic

K = 4;

Kx = K;
Ky = K;
Kz = K;

p = 4;

el_bounds.x = linspace(-1,1,Kx+1);
el_bounds.y = linspace(-1,1,Ky+1);
el_bounds.z = linspace(-1,1,Kz+1);

element.nodes = @(xii,eta,zta,elx,ely,elz) mesh.dim_3.element.mapping_nodes(xii,eta,zta,elx,ely,elz,el_bounds);

ttl_nr_el = Kx * Ky * Kz;

%% get discretization nodes

% dtize.nodes = get_dtize_nodes_3D();

%% get jacobian terms

dtize.jacobian = get_jacobian_3D(p,Kx,Ky,Kz,domain.jacobian); 

%% evaluate permeability atintegration nodes

perm.kxx = 1;
perm.kxy = 0;
perm.kxz = 0;
perm.kyy = 1;
perm.kyz = 0;
perm.kzz = 1;

%% get metric terms

dtize.metric = metric.dim_3.get_metric_3D_darcy(dtize.jacobian,perm);

%% create discretization nodes & basis

[x,w] = GLLnodes(p);
[~,e] = MimeticpolyVal(x,p,1);

h = eye(p+1); %%%%% this is only valid if you are using GLL Nodes - please check before using 

dtize.basis.hx = h;
dtize.basis.hy = h;
dtize.basis.hz = h;

dtize.basis.ex = e;
dtize.basis.ey = e;
dtize.basis.ez = e;

dtize.weights.wx = w;
dtize.weights.wy = w;
dtize.weights.wz = w;

dtize.nodes.sx = x;
dtize.nodes.sy = x;
dtize.nodes.sz = x;

dtize.metric.www = kron(kron(dtize.weights.wz,dtize.weights.wy),dtize.weights.wx);

%% create mass matrix 

dtize.basis.surface.yz = kron(kron(dtize.basis.ez,dtize.basis.ey),dtize.basis.hx);
dtize.basis.surface.zx = kron(kron(dtize.basis.ez,dtize.basis.hy),dtize.basis.ex);
dtize.basis.surface.xy = kron(kron(dtize.basis.hz,dtize.basis.ey),dtize.basis.ex);

M22 = get_stiffness_matrix_3D_surfaces(dtize.basis.surface,dtize.metric,ttl_nr_el);

%% random variables

nr_surfc = 3*p^2*(p+1);
nr_volum = p^3;
ttl_loc_dof = nr_surfc+nr_volum;
nr_local_lambda = 6*p^2;

nr_lambda = (Kx+1)*Ky*Kz*p^2 + (Ky+1)*Kx*Kz*p^2 + (Kz+1)*Kx*Ky*p^2;

%% create discretization nodes = reference

[dtize.nodes.eta, dtize.nodes.xii, dtize.nodes.zta] = meshgrid(dtize.nodes.sx,dtize.nodes.sy,dtize.nodes.sz);

dtize.nodes.xii = dtize.nodes.xii(:)';
dtize.nodes.eta = dtize.nodes.eta(:)';
dtize.nodes.zta = dtize.nodes.zta(:)';

%% create discretization nodes = physical

for i = 1:Kx
    for j = 1:Ky
        for k = 1:Kz
            eleid = (k-1)*Kx*Ky + (j-1)*Kx + i;
            [dtize.nodes.xii2(:,:,eleid),dtize.nodes.eta2(:,:,eleid),dtize.nodes.zta2(:,:,eleid)] = element.nodes(dtize.nodes.xii, dtize.nodes.eta, dtize.nodes.zta,i,j,k);
        end
    end
end

[dtize.nodes.x,dtize.nodes.y,dtize.nodes.z] = domain.mapping(dtize.nodes.xii2, dtize.nodes.eta2, dtize.nodes.zta2);

%% LHS

[gather_qqq] = gather_matrix.dim_3.continuous_surfaces(Kx, Ky, Kz, p);

for eln = 1:ttl_nr_el
    gather_ppp(1:nr_volum,eln) = (eln-1)*p^3 + (1:p^3);
end

E32 = full(incidence_matrix.dim_3.discrete_divergence(p));
E33 = repmat(E32, [1 1 ttl_nr_el]);

coll_M11 = AssembleMatrices.AssembleMatrices2(gather_qqq, gather_qqq, M22);
coll_E32 = AssembleMatrices.AssembleMatrices2(gather_ppp, gather_qqq, E33);

LHS = [coll_M11 coll_E32'; coll_E32 zeros(ttl_nr_el*p^3)];

% figure
% spy(LHS)

% xlabel('$$N$$', 'Interpreter','latex','FontSize',14);
% ylabel('$$\Vert \nabla \cdot \mathbf{q}^h - f^h \Vert _{L^2(\Omega)}$$', 'Interpreter','latex','FontSize',14);
% 
% set(gca,'TickLabelInterpreter','latex','FontSize',20,'Xcolor','k','Ycolor','k');

% export_fig('l2norm_divu_p_convergence.pdf','-pdf','-r864','-painters','-transparent');

% return

%% reduction of f

basis_volume = kron(kron(dtize.basis.ez,dtize.basis.ey),dtize.basis.ex);

dof_f = get_dof_3D_volumes(test.darcy_fff,basis_volume,dtize.metric, dtize.nodes, nr_volum, ttl_nr_el);

%% dirichlet boundary conditions

%% left face

basis_left_face = kron(dtize.basis.ez,dtize.basis.ey);

[nodes_left_face.xii,nodes_left_face.eta,nodes_left_face.zta] = meshgrid(-1,dtize.nodes.sy,dtize.nodes.sz);

nodes_left_face.xii = nodes_left_face.xii(:)';
nodes_left_face.eta = nodes_left_face.eta(:)';
nodes_left_face.zta = nodes_left_face.zta(:)';

% get new ref nodes
[nodes_left_face.xii2,nodes_left_face.eta2,nodes_left_face.zta2] = new_ref_nodes(1,1:Ky,1:Kz,nodes_left_face.xii,nodes_left_face.eta,nodes_left_face.zta,Ky,element.nodes);

[nodes_left_face.xxx,nodes_left_face.yyy,nodes_left_face.zzz] = domain.mapping(nodes_left_face.xii2,nodes_left_face.eta2,nodes_left_face.zta2);

phi_left_face = test.phi(nodes_left_face.xxx,nodes_left_face.yyy,nodes_left_face.zzz);

w88_left_face = kron(dtize.weights.wz,dtize.weights.wy);

nn_left = -1;

for eln = 1:Ky*Kz
    boundary_dof.left_face(:,eln) = nn_left * reduction(phi_left_face(:,:,eln), basis_left_face, w88_left_face, 1);
end

%% right face

basis_rght_face = kron(dtize.basis.ez,dtize.basis.ey);

[nodes_rght_face.xii,nodes_rght_face.eta,nodes_rght_face.zta] = meshgrid(+1,dtize.nodes.sy,dtize.nodes.sz);

nodes_rght_face.xii = nodes_rght_face.xii(:)';
nodes_rght_face.eta = nodes_rght_face.eta(:)';
nodes_rght_face.zta = nodes_rght_face.zta(:)';

% get new ref nodes
[nodes_rght_face.xii2,nodes_rght_face.eta2,nodes_rght_face.zta2] = new_ref_nodes(Kx,1:Ky,1:Kz,nodes_rght_face.xii,nodes_rght_face.eta,nodes_rght_face.zta,Ky,element.nodes);

[nodes_rght_face.xxx,nodes_rght_face.yyy,nodes_rght_face.zzz] = domain.mapping(nodes_rght_face.xii2,nodes_rght_face.eta2,nodes_rght_face.zta2);

phi_rght_face = test.phi(nodes_rght_face.xxx,nodes_rght_face.yyy,nodes_rght_face.zzz);

w88_rght_face = kron(dtize.weights.wz,dtize.weights.wy);

nn_rght = +1;

for eln = 1:Ky*Kz
    boundary_dof.rght_face(:,eln) = nn_rght * reduction(phi_rght_face(:,:,eln), basis_rght_face, w88_rght_face, 1); 
end

%% bottom face

basis_bttm_face = kron(dtize.basis.ez,dtize.basis.ex);
[nodes_bttm_face.xii,nodes_bttm_face.eta,nodes_bttm_face.zta] = meshgrid(dtize.nodes.sx,-1,dtize.nodes.sz);

nodes_bttm_face.xii = nodes_bttm_face.xii(:)';
nodes_bttm_face.eta = nodes_bttm_face.eta(:)';
nodes_bttm_face.zta = nodes_bttm_face.zta(:)';

% get new ref nodes
[nodes_bttm_face.xii2,nodes_bttm_face.eta2,nodes_bttm_face.zta2] = new_ref_nodes(1:Kx,1,1:Kz,nodes_bttm_face.xii,nodes_bttm_face.eta,nodes_bttm_face.zta,Kx,element.nodes);

for i = 1:Kx
    for k = 1:Kz
        local_eleid = (k-1)*Kx + i;
        [dtize.nodes_bttm_face.xii2(:,:,local_eleid),dtize.nodes_bttm_face.eta2(:,:,local_eleid),dtize.nodes_bttm_face.zta2(:,:,local_eleid)] = element.nodes(nodes_bttm_face.xii,nodes_bttm_face.eta,nodes_bttm_face.zta,i,1,k);
    end
end

w88_bttm_face = kron(dtize.weights.wz,dtize.weights.wx);

nn_bttm = -1;

boundary_dof.bttm_face = generate_boundary_dof(dtize.nodes_bttm_face, domain.mapping, test.phi, basis_bttm_face, w88_bttm_face, nn_bttm, Kx*Kz);

%% top face

basis_topp_face = kron(dtize.basis.ez,dtize.basis.ex);

[nodes_topp_face.xii,nodes_topp_face.eta,nodes_topp_face.zta] = meshgrid(dtize.nodes.sx,+1,dtize.nodes.sz);

nodes_topp_face.xii = nodes_topp_face.xii(:)';
nodes_topp_face.eta = nodes_topp_face.eta(:)';
nodes_topp_face.zta = nodes_topp_face.zta(:)';

for i = 1:Kx
    for k = 1:Kz
        local_eleid = (k-1)*Kx + i;
        [dtize.nodes_topp_face.xii2(:,:,local_eleid),dtize.nodes_topp_face.eta2(:,:,local_eleid),dtize.nodes_topp_face.zta2(:,:,local_eleid)] = element.nodes(nodes_topp_face.xii,nodes_topp_face.eta,nodes_topp_face.zta,i,Ky,k);
    end
end

w88_topp_face = kron(dtize.weights.wz,dtize.weights.wx);

nn_topp = +1;

boundary_dof.topp_face = generate_boundary_dof(dtize.nodes_topp_face, domain.mapping, test.phi, basis_topp_face, w88_topp_face, nn_topp, Kx*Kz);

%% back face

basis_back_face = kron(dtize.basis.ey,dtize.basis.ex);

[nodes_back_face.eta,nodes_back_face.xii,nodes_back_face.zta] = meshgrid(dtize.nodes.sx,dtize.nodes.sy,-1);

nodes_back_face.xii = nodes_back_face.xii(:)';
nodes_back_face.eta = nodes_back_face.eta(:)';
nodes_back_face.zta = nodes_back_face.zta(:)';

for i = 1:Kx
    for j = 1:Ky
        local_eleid = (j-1)*Kx + i;
        [nodes_back_face.xii2(:,:,local_eleid),nodes_back_face.eta2(:,:,local_eleid),nodes_back_face.zta2(:,:,local_eleid)] = element.nodes(nodes_back_face.xii,nodes_back_face.eta,nodes_back_face.zta,i,j,1);
    end
end

% [nodes_back_face.xxx,nodes_back_face.yyy,nodes_back_face.zzz] = domain.mapping(nodes_back_face.xii2,nodes_back_face.eta2,nodes_back_face.zta2);

% phi_back_face = test.phi(nodes_back_face.xxx,nodes_back_face.yyy,nodes_back_face.zzz);

w88_back_face = kron(dtize.weights.wy,dtize.weights.wx);

nn_back = -1;

boundary_dof.back_face = generate_boundary_dof(nodes_back_face, domain.mapping, test.phi, basis_back_face, w88_back_face, nn_back, Kx*Ky);

%% front face

basis_frnt_face = kron(dtize.basis.ey,dtize.basis.ex);

[nodes_frnt_face.eta,nodes_frnt_face.xii,nodes_frnt_face.zta] = meshgrid(dtize.nodes.sx,dtize.nodes.sy,+1);

nodes_frnt_face.xii = nodes_frnt_face.xii(:)';
nodes_frnt_face.eta = nodes_frnt_face.eta(:)';
nodes_frnt_face.zta = nodes_frnt_face.zta(:)';

for i = 1:Kx
    for j = 1:Ky
        local_eleid = (j-1)*Kx + i;
        [nodes_frnt_face.xii2(:,:,local_eleid),nodes_frnt_face.eta2(:,:,local_eleid),nodes_frnt_face.zta2(:,:,local_eleid)] = element.nodes(nodes_frnt_face.xii,nodes_frnt_face.eta,nodes_frnt_face.zta,i,j,Kz);
    end
end

% [nodes_frnt_face.xxx,nodes_frnt_face.yyy,nodes_frnt_face.zzz] = domain.mapping(nodes_frnt_face.xii2,nodes_frnt_face.eta2,nodes_frnt_face.zta2);

% phi_frnt_face = test.phi(nodes_frnt_face.xxx,nodes_frnt_face.yyy,nodes_frnt_face.zzz);

w88_frnt_face = kron(dtize.weights.wy,dtize.weights.wx);

nn_frnt = +1;

boundary_dof.frnt_face = generate_boundary_dof(nodes_frnt_face, domain.mapping, test.phi, basis_frnt_face, w88_frnt_face, nn_frnt, Kx*Ky);

%% add Dirichlet BCs' to RHS 2

element_face = mesh.dim_3.element_boundary_faces_map(p);

% left 
for j = 1:Ky
    for k = 1:Kz
        local_eleid = (k-1)*Ky + j;
        globl_eleid = (k-1)*Kx*Ky + (j-1)*Kx + 1;
        local_left_boundary_face = element_face.left(1:p^2);
        global_id = gather_qqq(local_left_boundary_face,globl_eleid);
        RHSX(global_id) = boundary_dof.left_face(:,local_eleid);
    end
end

% right
for j = 1:Ky
    for k = 1:Kz
        local_eleid = (k-1)*Ky + j;
        globl_eleid = (k-1)*Kx*Ky + (j-1)*Kx + Kx;
        local_rght_boundary_face = element_face.rght(1:p^2);
        global_id = gather_qqq(local_rght_boundary_face,globl_eleid);
        RHSX(global_id) = boundary_dof.rght_face(:,local_eleid);
    end
end

% bottom 
for i = 1:Kx
    for k = 1:Kz
        local_eleid = (k-1)*Kx + i;
        globl_eleid = (k-1)*Kx*Ky + (1-1)*Kx + i;
        local_bttm_boundary_face = element_face.bttm(1:p^2);
        global_id = gather_qqq(local_bttm_boundary_face,globl_eleid);
        RHSX(global_id) = boundary_dof.bttm_face(:,local_eleid);
    end
end

% top
for i = 1:Kx
    for k = 1:Kz
        local_eleid = (k-1)*Kx + i;
        globl_eleid = (k-1)*Kx*Ky + (Ky-1)*Kx + i;
        local_topp_boundary_face = element_face.topp(1:p^2);
        global_id = gather_qqq(local_topp_boundary_face,globl_eleid);
        RHSX(global_id) = boundary_dof.topp_face(:,local_eleid);
    end
end

% back
for i = 1:Kx
    for j = 1:Ky
        local_eleid = (j-1)*Kx + i;
        globl_eleid = (1-1)*Kx*Ky + (j-1)*Kx + i;
        local_back_boundary_face = element_face.back(1:p^2);
        global_id = gather_qqq(local_back_boundary_face,globl_eleid);
        RHSX(global_id) = boundary_dof.back_face(:,local_eleid);
    end
end

% front
for i = 1:Kx
    for j = 1:Ky
        local_eleid = (j-1)*Kx + i;
        globl_eleid = (Kz-1)*Kx*Ky + (j-1)*Kx + i;
        local_frnt_boundary_face = element_face.frnt(1:p^2);
        global_id = gather_qqq(local_frnt_boundary_face,globl_eleid);
        RHSX(global_id) = boundary_dof.frnt_face(:,local_eleid);
    end
end

%% RHS

RHSp = dof_f(:);

RHSn = [RHSX'; RHSp];

%% solution 

X = LHS\RHSn;

% X = cgs(LHS,RHSn,1e-20,1000);

toc

% return

%% reconstruction

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

%% create reconstruction jacobian and nodes

[rec.nodes.eta,rec.nodes.xii,rec.nodes.zta] = meshgrid(rec.nodes.sfx,rec.nodes.sfy,rec.nodes.sfz);

rec.nodes.xii = rec.nodes.xii(:)';
rec.nodes.eta = rec.nodes.eta(:)';
rec.nodes.zta = rec.nodes.zta(:)';

for i = 1:Kx
    for j = 1:Ky
        for k = 1:Kz
            eleid = (k-1)*Kx*Ky + (j-1)*Kx + i;
            [rec.nodes.xii2(:,:,eleid),rec.nodes.eta2(:,:,eleid),rec.nodes.zta2(:,:,eleid)] = element.nodes(rec.nodes.xii, rec.nodes.eta, rec.nodes.zta,i,j,k);
        end
    end
end

[rec.nodes.x,rec.nodes.y,rec.nodes.z] = domain.mapping(rec.nodes.xii2,rec.nodes.eta2,rec.nodes.zta2);

rec.jacobian = domain.jacobian(rec.nodes.xii2, rec.nodes.eta2, rec.nodes.zta2);

dXii_ds = 0.5*(el_bounds.x(2:end) - el_bounds.x(1:end-1));
dEta_dt = 0.5*(el_bounds.y(2:end) - el_bounds.y(1:end-1));
dZta_du = 0.5*(el_bounds.z(2:end) - el_bounds.z(1:end-1));

[d1, d2, d3] = meshgrid(dXii_ds,dEta_dt,dZta_du);

d1 = permute(d1,[2 1 3]);
d1 = d1(:)';

d2 = permute(d2,[2 1 3]);
d2 = d2(:)';

d3 = d3(:)';

rec.jacobian.dXdxii = rec.jacobian.dXdxii .* reshape(d1, [1 1 ttl_nr_el]);
rec.jacobian.dYdxii = rec.jacobian.dYdxii .* reshape(d1, [1 1 ttl_nr_el]);
rec.jacobian.dZdxii = rec.jacobian.dZdxii .* reshape(d1, [1 1 ttl_nr_el]);

rec.jacobian.dXdeta = rec.jacobian.dXdeta .* reshape(d2, [1 1 ttl_nr_el]);
rec.jacobian.dYdeta = rec.jacobian.dYdeta .* reshape(d2, [1 1 ttl_nr_el]);
rec.jacobian.dZdeta = rec.jacobian.dZdeta .* reshape(d2, [1 1 ttl_nr_el]);

rec.jacobian.dXdzta = rec.jacobian.dXdzta .* reshape(d3, [1 1 ttl_nr_el]);
rec.jacobian.dYdzta = rec.jacobian.dYdzta .* reshape(d3, [1 1 ttl_nr_el]);
rec.jacobian.dZdzta = rec.jacobian.dZdzta .* reshape(d3, [1 1 ttl_nr_el]);

% for eln = 1:ttl_nr_el
%     rec.jacobian.dXdxii(:,:,eln) = rec.jacobian.dXdxii(:,:,eln) * d1(eln);
%     rec.jacobian.dYdxii(:,:,eln) = rec.jacobian.dYdxii(:,:,eln) * d1(eln);
%     rec.jacobian.dZdxii(:,:,eln) = rec.jacobian.dZdxii(:,:,eln) * d1(eln);
%     
%     rec.jacobian.dXdeta(:,:,eln) = rec.jacobian.dXdeta(:,:,eln) * d2(eln);
%     rec.jacobian.dYdeta(:,:,eln) = rec.jacobian.dYdeta(:,:,eln) * d2(eln);
%     rec.jacobian.dZdeta(:,:,eln) = rec.jacobian.dZdeta(:,:,eln) * d2(eln);
%     
%     rec.jacobian.dXdzta(:,:,eln) = rec.jacobian.dXdzta(:,:,eln) * d3(eln);
%     rec.jacobian.dYdzta(:,:,eln) = rec.jacobian.dYdzta(:,:,eln) * d3(eln);
%     rec.jacobian.dZdzta(:,:,eln) = rec.jacobian.dZdzta(:,:,eln) * d3(eln);
% end

% determinant of jacobian
rec.metric.ggg = metric.dim_3.determinant_jacobian(rec.jacobian);
rec.metric.invg = 1 ./rec.metric.ggg;

rec.wwwf = kron(kron(rec.weights.wfz,rec.weights.wfy),rec.weights.wfx);

%% create reconstruction basis and nodes - II

rec_basis_volume = kron(kron(rec.basis.efz,rec.basis.efy),rec.basis.efx);

%% reconstruction of f

rec.metric.ggg = metric.dim_3.determinant_jacobian(rec.jacobian);
rec.metric.invg = 1 ./rec.metric.ggg;

for eln = 1:ttl_nr_el
    rec.fff(:,:,eln) = reconstruction.dim_3.volume(dof_f(:,1,eln),rec_basis_volume,pf,rec.metric.ggg(:,:,eln))';
end

ext_fff = test.darcy_fff(rec.nodes.x,rec.nodes.y,rec.nodes.z);

error_fff = sqrt(sum(sum((ext_fff-rec.fff).^2 .* rec.metric.ggg .* rec.wwwf)));

error_fff

%% reduction of q^(n-1)

qqx_int = test.dpdx(dtize.nodes.x,dtize.nodes.y,dtize.nodes.z);
qqy_int = test.dpdy(dtize.nodes.x,dtize.nodes.y,dtize.nodes.z);
qqz_int = test.dpdz(dtize.nodes.x,dtize.nodes.y,dtize.nodes.z);

%% reconstruction of q^(n-1)

rec.basis.surfac.xy = kron(kron(rec.basis.hfz,rec.basis.efy),rec.basis.efx);
rec.basis.surfac.yz = kron(kron(rec.basis.efz,rec.basis.efy),rec.basis.hfx);
rec.basis.surfac.zx = kron(kron(rec.basis.efz,rec.basis.hfy),rec.basis.efx);

%% cochains

X_u = X(1:3*(K*p)^2*((K*p)+1));

dof_u = X_u(gather_qqq);

X_p = X(3*(K*p)^2*((K*p)+1)+1:3*(K*p)^2*((K*p)+1)+(K*p)^3);

dof_p = X_p(gather_ppp);

if p==1
    dof_p = dof_p';
end

% dof_p = X(nr_surfc+1:nr_surfc+nr_volum,:);

%% reconstruction of div q^h 

for eln = 1:ttl_nr_el
    divQ_h(:,eln) = E32 * dof_u(:,eln);
end

for eln = 1:ttl_nr_el
    rec.divQ(:,:,eln) = reconstruction.dim_3.volume(divQ_h(:,eln),rec_basis_volume,pf,rec.metric.ggg(:,:,eln))';
end

error_divQ_h = sqrt(sum(sum((ext_fff-rec.divQ).^2 .* rec.metric.ggg .* rec.wwwf)));

error_divQ_h

%% reconstruction of constraint div q^h - f^h 

for eln = 1:ttl_nr_el
    dof_f2(:,eln) = dof_f(:,1,eln);
end

err_const_h = divQ_h - dof_f2;

for eln = 1:ttl_nr_el
    rec.const(:,:,eln) = reconstruction.dim_3.volume(err_const_h(:,eln),rec_basis_volume,pf,rec.metric.ggg(:,:,eln))';
end

error_const_h = sqrt(sum(sum((rec.const).^2 .* rec.metric.ggg .* rec.wwwf)));

error_const_h

%% 

M333 = zeros(nr_volum, nr_volum, ttl_nr_el);

for eln = 1:ttl_nr_el
    M333(:,:,eln) = construct_mass_matrix(basis_volume,basis_volume,dtize.metric.www,dtize.metric.invg(:,:,eln));
end


%%

for eln = 1:ttl_nr_el
    dof_p(:,eln) = M333(:,:,eln) \ dof_p(:,eln);
end

%% reconstruction pressure

for eln = 1:ttl_nr_el
    rec.ppp(:,:,eln) = reconstruction.dim_3.volume(dof_p(:,eln),rec_basis_volume,pf,rec.metric.ggg(:,:,eln))';
end

ext_ppp = test.phi(rec.nodes.x,rec.nodes.y,rec.nodes.z);
ext_qxx = test.darcy_qx(rec.nodes.x,rec.nodes.y,rec.nodes.z);
ext_qyy = test.darcy_qy(rec.nodes.x,rec.nodes.y,rec.nodes.z);
ext_qzz = test.darcy_qz(rec.nodes.x,rec.nodes.y,rec.nodes.z);

error_ppp = sqrt(sum(sum((ext_ppp-rec.ppp).^2 .* rec.metric.ggg .* rec.wwwf)));

% error_fff
error_ppp

% return

%% plot the pressure

element_nr = @(elx,ely,elz) (elz-1)*Kx*Ky + (ely-1)*Kx + elx;

% find boundary elements 

% back/frnt boundary 
back_eleid = 1:Kx*Ky;
frnt_eleid = (Kz-1)*Kx*Ky + back_eleid;

%left / right
for j = 1:Ky
    for k =1:Kz
        temp = (k-1)*Ky + j;
        left_eleid(temp) = element_nr(1,j,k);
        rght_eleid(temp) = element_nr(Kx,j,k);
    end
end

% bttm / top

for i = 1:Kx
    for k = 1:Kz
        temp = (k-1)*Kx + i;
        bttm_eleid(temp) = element_nr(i,1,k);
        topp_eleid(temp) = element_nr(i,Ky,k);
    end
end

bndry_eleid = [back_eleid frnt_eleid left_eleid rght_eleid bttm_eleid topp_eleid];
ttl_nr_boundary_faces = 2*Kx*Ky + 2*Ky*Kz + 2*Kx*Kz;

figure
hold on 

% lift right face
for eln2 = 1:Ky*Kz
    sliceid = 1;
    eln = left_eleid(eln2);
    
    plt_ext_ppp = reshape(ext_ppp(:,:,eln), [pf+1 pf+1 pf+1]);
    plt_rec_xxx = reshape(rec.nodes.x(:,:,eln), [pf+1 pf+1 pf+1]);
    plt_rec_yyy = reshape(rec.nodes.y(:,:,eln), [pf+1 pf+1 pf+1]);
    plt_rec_zzz = reshape(rec.nodes.z(:,:,eln), [pf+1 pf+1 pf+1]);

    surf(reshape(plt_rec_xxx(sliceid,:,:),[pf+1,pf+1]),reshape(plt_rec_yyy(sliceid,:,:),[pf+1,pf+1]),reshape(plt_rec_zzz(sliceid,:,:),[pf+1,pf+1]),reshape(plt_ext_ppp(sliceid,:,:),[pf+1,pf+1]));
%     set(h, 'EdgeColor','none')
    
    sliceid = pf+1;
    eln = rght_eleid(eln2);
    
    plt_ext_ppp = reshape(ext_ppp(:,:,eln), [pf+1 pf+1 pf+1]);
    plt_rec_xxx = reshape(rec.nodes.x(:,:,eln), [pf+1 pf+1 pf+1]);
    plt_rec_yyy = reshape(rec.nodes.y(:,:,eln), [pf+1 pf+1 pf+1]);
    plt_rec_zzz = reshape(rec.nodes.z(:,:,eln), [pf+1 pf+1 pf+1]);

    surf(reshape(plt_rec_xxx(sliceid,:,:),[pf+1,pf+1]),reshape(plt_rec_yyy(sliceid,:,:),[pf+1,pf+1]),reshape(plt_rec_zzz(sliceid,:,:),[pf+1,pf+1]),reshape(plt_ext_ppp(sliceid,:,:),[pf+1,pf+1]));
%     set(h, 'EdgeColor','none')
end

% top bottom face 
for eln2 = 1:Kx*Kz
    sliceid = 1;
    eln = bttm_eleid(eln2);
    
    plt_ext_ppp = reshape(ext_ppp(:,:,eln), [pf+1 pf+1 pf+1]);
    plt_rec_xxx = reshape(rec.nodes.x(:,:,eln), [pf+1 pf+1 pf+1]);
    plt_rec_yyy = reshape(rec.nodes.y(:,:,eln), [pf+1 pf+1 pf+1]);
    plt_rec_zzz = reshape(rec.nodes.z(:,:,eln), [pf+1 pf+1 pf+1]);

    surf(reshape(plt_rec_xxx(:,sliceid,:),[pf+1,pf+1]),reshape(plt_rec_yyy(:,sliceid,:),[pf+1,pf+1]),reshape(plt_rec_zzz(:,sliceid,:),[pf+1,pf+1]),reshape(plt_ext_ppp(:,sliceid,:),[pf+1,pf+1]));
%     set(h, 'EdgeColor','none')
    
    sliceid = pf+1;
    eln = topp_eleid(eln2);
    
    plt_ext_ppp = reshape(ext_ppp(:,:,eln), [pf+1 pf+1 pf+1]);
    plt_rec_xxx = reshape(rec.nodes.x(:,:,eln), [pf+1 pf+1 pf+1]);
    plt_rec_yyy = reshape(rec.nodes.y(:,:,eln), [pf+1 pf+1 pf+1]);
    plt_rec_zzz = reshape(rec.nodes.z(:,:,eln), [pf+1 pf+1 pf+1]);

    surf(reshape(plt_rec_xxx(:,sliceid,:),[pf+1,pf+1]),reshape(plt_rec_yyy(:,sliceid,:),[pf+1,pf+1]),reshape(plt_rec_zzz(:,sliceid,:),[pf+1,pf+1]),reshape(plt_ext_ppp(:,sliceid,:),[pf+1,pf+1]));
%     set(h, 'EdgeColor','none')
end

% front back face
for eln2 = 1:Kx*Ky
    sliceid = 1;
    eln = back_eleid(eln2);
    
    plt_ext_ppp = reshape(ext_ppp(:,:,eln), [pf+1 pf+1 pf+1]);
    plt_rec_xxx = reshape(rec.nodes.x(:,:,eln), [pf+1 pf+1 pf+1]);
    plt_rec_yyy = reshape(rec.nodes.y(:,:,eln), [pf+1 pf+1 pf+1]);
    plt_rec_zzz = reshape(rec.nodes.z(:,:,eln), [pf+1 pf+1 pf+1]);

    surf(reshape(plt_rec_xxx(:,:,sliceid),[pf+1,pf+1]),reshape(plt_rec_yyy(:,:,sliceid),[pf+1,pf+1]),reshape(plt_rec_zzz(:,:,sliceid),[pf+1,pf+1]),reshape(plt_ext_ppp(:,:,sliceid),[pf+1,pf+1]));
%     set(h, 'EdgeColor','none')
    
    sliceid = pf+1;
    eln = frnt_eleid(eln2);
    
    plt_ext_ppp = reshape(ext_ppp(:,:,eln), [pf+1 pf+1 pf+1]);
    plt_rec_xxx = reshape(rec.nodes.x(:,:,eln), [pf+1 pf+1 pf+1]);
    plt_rec_yyy = reshape(rec.nodes.y(:,:,eln), [pf+1 pf+1 pf+1]);
    plt_rec_zzz = reshape(rec.nodes.z(:,:,eln), [pf+1 pf+1 pf+1]);

    surf(reshape(plt_rec_xxx(:,:,sliceid),[pf+1,pf+1]),reshape(plt_rec_yyy(:,:,sliceid),[pf+1,pf+1]),reshape(plt_rec_zzz(:,:,sliceid),[pf+1,pf+1]),reshape(plt_ext_ppp(:,:,sliceid),[pf+1,pf+1]));
%     set(h, 'EdgeColor','none')
end

colorbar
title('pressure exact')
colormap jet

figure
hold on 

% lift right face
for eln2 = 1:Ky*Kz
    sliceid = 1;
    eln = left_eleid(eln2);
    
    plt_rec_ppp = reshape(rec.ppp(:,:,eln), [pf+1 pf+1 pf+1]);
    plt_rec_xxx = reshape(rec.nodes.x(:,:,eln), [pf+1 pf+1 pf+1]);
    plt_rec_yyy = reshape(rec.nodes.y(:,:,eln), [pf+1 pf+1 pf+1]);
    plt_rec_zzz = reshape(rec.nodes.z(:,:,eln), [pf+1 pf+1 pf+1]);

    surf(reshape(plt_rec_xxx(sliceid,:,:),[pf+1,pf+1]),reshape(plt_rec_yyy(sliceid,:,:),[pf+1,pf+1]),reshape(plt_rec_zzz(sliceid,:,:),[pf+1,pf+1]),reshape(plt_rec_ppp(sliceid,:,:),[pf+1,pf+1]));
%     set(h, 'EdgeColor','none')
    
    sliceid = pf+1;
    eln = rght_eleid(eln2);
    
    plt_rec_ppp = reshape(rec.ppp(:,:,eln), [pf+1 pf+1 pf+1]);
    plt_rec_xxx = reshape(rec.nodes.x(:,:,eln), [pf+1 pf+1 pf+1]);
    plt_rec_yyy = reshape(rec.nodes.y(:,:,eln), [pf+1 pf+1 pf+1]);
    plt_rec_zzz = reshape(rec.nodes.z(:,:,eln), [pf+1 pf+1 pf+1]);

    surf(reshape(plt_rec_xxx(sliceid,:,:),[pf+1,pf+1]),reshape(plt_rec_yyy(sliceid,:,:),[pf+1,pf+1]),reshape(plt_rec_zzz(sliceid,:,:),[pf+1,pf+1]),reshape(plt_rec_ppp(sliceid,:,:),[pf+1,pf+1]));
%     set(h, 'EdgeColor','none')
end

% top bottom face 
for eln2 = 1:Kx*Kz
    sliceid = 1;
    eln = bttm_eleid(eln2);
    
    plt_rec_ppp = reshape(rec.ppp(:,:,eln), [pf+1 pf+1 pf+1]);
    plt_rec_xxx = reshape(rec.nodes.x(:,:,eln), [pf+1 pf+1 pf+1]);
    plt_rec_yyy = reshape(rec.nodes.y(:,:,eln), [pf+1 pf+1 pf+1]);
    plt_rec_zzz = reshape(rec.nodes.z(:,:,eln), [pf+1 pf+1 pf+1]);

    surf(reshape(plt_rec_xxx(:,sliceid,:),[pf+1,pf+1]),reshape(plt_rec_yyy(:,sliceid,:),[pf+1,pf+1]),reshape(plt_rec_zzz(:,sliceid,:),[pf+1,pf+1]),reshape(plt_rec_ppp(:,sliceid,:),[pf+1,pf+1]));
%     set(h, 'EdgeColor','none')
    
    sliceid = pf+1;
    eln = topp_eleid(eln2);
    
    plt_rec_ppp = reshape(rec.ppp(:,:,eln), [pf+1 pf+1 pf+1]);
    plt_rec_xxx = reshape(rec.nodes.x(:,:,eln), [pf+1 pf+1 pf+1]);
    plt_rec_yyy = reshape(rec.nodes.y(:,:,eln), [pf+1 pf+1 pf+1]);
    plt_rec_zzz = reshape(rec.nodes.z(:,:,eln), [pf+1 pf+1 pf+1]);

    surf(reshape(plt_rec_xxx(:,sliceid,:),[pf+1,pf+1]),reshape(plt_rec_yyy(:,sliceid,:),[pf+1,pf+1]),reshape(plt_rec_zzz(:,sliceid,:),[pf+1,pf+1]),reshape(plt_rec_ppp(:,sliceid,:),[pf+1,pf+1]));
%     set(h, 'EdgeColor','none')
end

% front back face
for eln2 = 1:Kx*Ky
    sliceid = 1;
    eln = back_eleid(eln2);
    
    plt_rec_ppp = reshape(rec.ppp(:,:,eln), [pf+1 pf+1 pf+1]);
    plt_rec_xxx = reshape(rec.nodes.x(:,:,eln), [pf+1 pf+1 pf+1]);
    plt_rec_yyy = reshape(rec.nodes.y(:,:,eln), [pf+1 pf+1 pf+1]);
    plt_rec_zzz = reshape(rec.nodes.z(:,:,eln), [pf+1 pf+1 pf+1]);

    surf(reshape(plt_rec_xxx(:,:,sliceid),[pf+1,pf+1]),reshape(plt_rec_yyy(:,:,sliceid),[pf+1,pf+1]),reshape(plt_rec_zzz(:,:,sliceid),[pf+1,pf+1]),reshape(plt_rec_ppp(:,:,sliceid),[pf+1,pf+1]));
%     set(h, 'EdgeColor','none')
    
    sliceid = pf+1;
    eln = frnt_eleid(eln2);
    
    plt_rec_ppp = reshape(rec.ppp(:,:,eln), [pf+1 pf+1 pf+1]);
    plt_rec_xxx = reshape(rec.nodes.x(:,:,eln), [pf+1 pf+1 pf+1]);
    plt_rec_yyy = reshape(rec.nodes.y(:,:,eln), [pf+1 pf+1 pf+1]);
    plt_rec_zzz = reshape(rec.nodes.z(:,:,eln), [pf+1 pf+1 pf+1]);

    surf(reshape(plt_rec_xxx(:,:,sliceid),[pf+1,pf+1]),reshape(plt_rec_yyy(:,:,sliceid),[pf+1,pf+1]),reshape(plt_rec_zzz(:,:,sliceid),[pf+1,pf+1]),reshape(plt_rec_ppp(:,:,sliceid),[pf+1,pf+1]));
%     set(h, 'EdgeColor','none')
end

colorbar
title('pressure reconstruction')
colormap jet

%% return 

% return

%% reconstruction flux

figure
hold on 

for eln = 1:ttl_nr_el
    plt_ext_qxx = reshape(ext_qxx(:,:,eln), [pf+1 pf+1 pf+1]);
    plt_rec_xxx = reshape(rec.nodes.x(:,:,eln), [pf+1 pf+1 pf+1]);
    plt_rec_yyy = reshape(rec.nodes.y(:,:,eln), [pf+1 pf+1 pf+1]);
    plt_rec_zzz = reshape(rec.nodes.z(:,:,eln), [pf+1 pf+1 pf+1]);

    surf(reshape(plt_rec_xxx(sliceid,:,:),[pf+1,pf+1]),reshape(plt_rec_yyy(sliceid,:,:),[pf+1,pf+1]),reshape(plt_rec_zzz(sliceid,:,:),[pf+1,pf+1]),reshape(plt_ext_qxx(sliceid,:,:),[pf+1,pf+1]))
    surf(reshape(plt_rec_xxx(:,sliceid,:),[pf+1,pf+1]),reshape(plt_rec_yyy(:,sliceid,:),[pf+1,pf+1]),reshape(plt_rec_zzz(:,sliceid,:),[pf+1,pf+1]),reshape(plt_ext_qxx(:,sliceid,:),[pf+1,pf+1]))
    surf(reshape(plt_rec_xxx(:,:,sliceid),[pf+1,pf+1]),reshape(plt_rec_yyy(:,:,sliceid),[pf+1,pf+1]),reshape(plt_rec_zzz(:,:,sliceid),[pf+1,pf+1]),reshape(plt_ext_qxx(:,:,sliceid),[pf+1,pf+1]))
end

colorbar
title('velocity x exact')

for eln = 1:ttl_nr_el
%     [rec.u.ux(:,:,eln), rec.u.uy(:,:,eln), rec.u.uz(:,:,eln)] = reconstruction.dim_3.surface(X(1:nr_surfc,eln),rec.jacobian,p, rec.metric.invg, rec.basis.surfac);
    
    local_jac.dXdxii = rec.jacobian.dXdxii(:,:,eln);
    local_jac.dXdeta = rec.jacobian.dXdeta(:,:,eln);
    local_jac.dXdzta = rec.jacobian.dXdzta(:,:,eln);
    local_jac.dYdxii = rec.jacobian.dYdxii(:,:,eln);
    local_jac.dYdeta = rec.jacobian.dYdeta(:,:,eln);
    local_jac.dYdzta = rec.jacobian.dYdzta(:,:,eln);
    local_jac.dZdxii = rec.jacobian.dZdxii(:,:,eln);
    local_jac.dZdeta = rec.jacobian.dZdeta(:,:,eln);
    local_jac.dZdzta = rec.jacobian.dZdzta(:,:,eln);
             
    [rec.u.ux(:,:,eln), rec.u.uy(:,:,eln), rec.u.uz(:,:,eln)] = reconstruction.dim_3.surface2(dof_u(:,eln),local_jac,p,rec.metric.invg(:,:,eln),rec.basis.surfac);
end

%% error flux

error_qxx = sqrt(sum(sum((ext_qxx-rec.u.ux).^2 .* rec.metric.ggg .* rec.wwwf)));
error_qyy = sqrt(sum(sum((ext_qyy-rec.u.uy).^2 .* rec.metric.ggg .* rec.wwwf)));
error_qzz = sqrt(sum(sum((ext_qzz-rec.u.uz).^2 .* rec.metric.ggg .* rec.wwwf)));

error_qxx
error_qyy
error_qzz

error_u_Hdiv = sqrt(sum(sum((ext_qxx-rec.u.ux).^2 .* rec.metric.ggg .* rec.wwwf)) ...
             + sum(sum((ext_qyy-rec.u.uy).^2 .* rec.metric.ggg .* rec.wwwf)) ...
             + sum(sum((ext_qzz-rec.u.uz).^2 .* rec.metric.ggg .* rec.wwwf)) ...
             + sum(sum(sum((rec.const).^2 .* rec.metric.ggg .* rec.wwwf))));

error_u_Hdiv

return

%% plots flux

figure
hold on

for eln = 1:ttl_nr_el
    plt_vel_xxx = reshape(rec.u.ux(:,:,eln), [pf+1 pf+1 pf+1]);
    plt_rec_xxx = reshape(rec.nodes.x(:,:,eln), [pf+1 pf+1 pf+1]);
    plt_rec_yyy = reshape(rec.nodes.y(:,:,eln), [pf+1 pf+1 pf+1]);
    plt_rec_zzz = reshape(rec.nodes.z(:,:,eln), [pf+1 pf+1 pf+1]);
    
    surf(reshape(plt_rec_xxx(sliceid,:,:),[pf+1,pf+1]),reshape(plt_rec_yyy(sliceid,:,:),[pf+1,pf+1]),reshape(plt_rec_zzz(sliceid,:,:),[pf+1,pf+1]),reshape(plt_vel_xxx(sliceid,:,:),[pf+1,pf+1]))
    surf(reshape(plt_rec_xxx(:,sliceid,:),[pf+1,pf+1]),reshape(plt_rec_yyy(:,sliceid,:),[pf+1,pf+1]),reshape(plt_rec_zzz(:,sliceid,:),[pf+1,pf+1]),reshape(plt_vel_xxx(:,sliceid,:),[pf+1,pf+1]))
    surf(reshape(plt_rec_xxx(:,:,sliceid),[pf+1,pf+1]),reshape(plt_rec_yyy(:,:,sliceid),[pf+1,pf+1]),reshape(plt_rec_zzz(:,:,sliceid),[pf+1,pf+1]),reshape(plt_vel_xxx(:,:,sliceid),[pf+1,pf+1]))
end

colorbar
title('velocity x reconstruction')

figure
hold on 

for eln = 1:ttl_nr_el
    plt_ext_qyy = reshape(ext_qyy(:,:,eln), [pf+1 pf+1 pf+1]);
    plt_rec_xxx = reshape(rec.nodes.x(:,:,eln), [pf+1 pf+1 pf+1]);
    plt_rec_yyy = reshape(rec.nodes.y(:,:,eln), [pf+1 pf+1 pf+1]);
    plt_rec_zzz = reshape(rec.nodes.z(:,:,eln), [pf+1 pf+1 pf+1]);

    surf(reshape(plt_rec_xxx(sliceid,:,:),[pf+1,pf+1]),reshape(plt_rec_yyy(sliceid,:,:),[pf+1,pf+1]),reshape(plt_rec_zzz(sliceid,:,:),[pf+1,pf+1]),reshape(plt_ext_qyy(sliceid,:,:),[pf+1,pf+1]))
    surf(reshape(plt_rec_xxx(:,sliceid,:),[pf+1,pf+1]),reshape(plt_rec_yyy(:,sliceid,:),[pf+1,pf+1]),reshape(plt_rec_zzz(:,sliceid,:),[pf+1,pf+1]),reshape(plt_ext_qyy(:,sliceid,:),[pf+1,pf+1]))
    surf(reshape(plt_rec_xxx(:,:,sliceid),[pf+1,pf+1]),reshape(plt_rec_yyy(:,:,sliceid),[pf+1,pf+1]),reshape(plt_rec_zzz(:,:,sliceid),[pf+1,pf+1]),reshape(plt_ext_qyy(:,:,sliceid),[pf+1,pf+1]))
end

colorbar
title('velocity y exact')

figure
hold on

for eln = 1:ttl_nr_el
    plt_vel_yyy = reshape(rec.u.uy(:,:,eln), [pf+1 pf+1 pf+1]);
    plt_rec_xxx = reshape(rec.nodes.x(:,:,eln), [pf+1 pf+1 pf+1]);
    plt_rec_yyy = reshape(rec.nodes.y(:,:,eln), [pf+1 pf+1 pf+1]);
    plt_rec_zzz = reshape(rec.nodes.z(:,:,eln), [pf+1 pf+1 pf+1]);
    
    surf(reshape(plt_rec_xxx(sliceid,:,:),[pf+1,pf+1]),reshape(plt_rec_yyy(sliceid,:,:),[pf+1,pf+1]),reshape(plt_rec_zzz(sliceid,:,:),[pf+1,pf+1]),reshape(plt_vel_yyy(sliceid,:,:),[pf+1,pf+1]))
    surf(reshape(plt_rec_xxx(:,sliceid,:),[pf+1,pf+1]),reshape(plt_rec_yyy(:,sliceid,:),[pf+1,pf+1]),reshape(plt_rec_zzz(:,sliceid,:),[pf+1,pf+1]),reshape(plt_vel_yyy(:,sliceid,:),[pf+1,pf+1]))
    surf(reshape(plt_rec_xxx(:,:,sliceid),[pf+1,pf+1]),reshape(plt_rec_yyy(:,:,sliceid),[pf+1,pf+1]),reshape(plt_rec_zzz(:,:,sliceid),[pf+1,pf+1]),reshape(plt_vel_yyy(:,:,sliceid),[pf+1,pf+1]))
end

colorbar
title('velocity y reconstruction')

figure
hold on 

for eln = 1:ttl_nr_el
    plt_ext_qzz = reshape(ext_qzz(:,:,eln), [pf+1 pf+1 pf+1]);
    plt_rec_xxx = reshape(rec.nodes.x(:,:,eln), [pf+1 pf+1 pf+1]);
    plt_rec_yyy = reshape(rec.nodes.y(:,:,eln), [pf+1 pf+1 pf+1]);
    plt_rec_zzz = reshape(rec.nodes.z(:,:,eln), [pf+1 pf+1 pf+1]);

    surf(reshape(plt_rec_xxx(sliceid,:,:),[pf+1,pf+1]),reshape(plt_rec_yyy(sliceid,:,:),[pf+1,pf+1]),reshape(plt_rec_zzz(sliceid,:,:),[pf+1,pf+1]),reshape(plt_ext_qzz(sliceid,:,:),[pf+1,pf+1]))
    surf(reshape(plt_rec_xxx(:,sliceid,:),[pf+1,pf+1]),reshape(plt_rec_yyy(:,sliceid,:),[pf+1,pf+1]),reshape(plt_rec_zzz(:,sliceid,:),[pf+1,pf+1]),reshape(plt_ext_qzz(:,sliceid,:),[pf+1,pf+1]))
    surf(reshape(plt_rec_xxx(:,:,sliceid),[pf+1,pf+1]),reshape(plt_rec_yyy(:,:,sliceid),[pf+1,pf+1]),reshape(plt_rec_zzz(:,:,sliceid),[pf+1,pf+1]),reshape(plt_ext_qzz(:,:,sliceid),[pf+1,pf+1]))
end

colorbar
title('velocity z exact')

figure
hold on

for eln = 1:ttl_nr_el
    plt_vel_zzz = reshape(rec.u.uz(:,:,eln), [pf+1 pf+1 pf+1]);
    plt_rec_xxx = reshape(rec.nodes.x(:,:,eln), [pf+1 pf+1 pf+1]);
    plt_rec_yyy = reshape(rec.nodes.y(:,:,eln), [pf+1 pf+1 pf+1]);
    plt_rec_zzz = reshape(rec.nodes.z(:,:,eln), [pf+1 pf+1 pf+1]);
    
    surf(reshape(plt_rec_xxx(sliceid,:,:),[pf+1,pf+1]),reshape(plt_rec_yyy(sliceid,:,:),[pf+1,pf+1]),reshape(plt_rec_zzz(sliceid,:,:),[pf+1,pf+1]),reshape(plt_vel_zzz(sliceid,:,:),[pf+1,pf+1]))
    surf(reshape(plt_rec_xxx(:,sliceid,:),[pf+1,pf+1]),reshape(plt_rec_yyy(:,sliceid,:),[pf+1,pf+1]),reshape(plt_rec_zzz(:,sliceid,:),[pf+1,pf+1]),reshape(plt_vel_zzz(:,sliceid,:),[pf+1,pf+1]))
    surf(reshape(plt_rec_xxx(:,:,sliceid),[pf+1,pf+1]),reshape(plt_rec_yyy(:,:,sliceid),[pf+1,pf+1]),reshape(plt_rec_zzz(:,:,sliceid),[pf+1,pf+1]),reshape(plt_vel_zzz(:,:,sliceid),[pf+1,pf+1]))
end

colorbar
title('velocity z reconstruction')

%%

figure
hold on

for eln = 1:ttl_nr_el
    plt_dif_qxx = reshape(ext_qxx(:,:,eln) - rec.u.ux(:,:,eln), [pf+1 pf+1 pf+1]);
    plt_rec_xxx = reshape(rec.nodes.x(:,:,eln), [pf+1 pf+1 pf+1]);
    plt_rec_yyy = reshape(rec.nodes.y(:,:,eln), [pf+1 pf+1 pf+1]);
    plt_rec_zzz = reshape(rec.nodes.z(:,:,eln), [pf+1 pf+1 pf+1]);
    
    surf(reshape(plt_rec_xxx(sliceid,:,:),[pf+1,pf+1]),reshape(plt_rec_yyy(sliceid,:,:),[pf+1,pf+1]),reshape(plt_rec_zzz(sliceid,:,:),[pf+1,pf+1]),reshape(plt_dif_qxx(sliceid,:,:),[pf+1,pf+1]))
    surf(reshape(plt_rec_xxx(:,sliceid,:),[pf+1,pf+1]),reshape(plt_rec_yyy(:,sliceid,:),[pf+1,pf+1]),reshape(plt_rec_zzz(:,sliceid,:),[pf+1,pf+1]),reshape(plt_dif_qxx(:,sliceid,:),[pf+1,pf+1]))
    surf(reshape(plt_rec_xxx(:,:,sliceid),[pf+1,pf+1]),reshape(plt_rec_yyy(:,:,sliceid),[pf+1,pf+1]),reshape(plt_rec_zzz(:,:,sliceid),[pf+1,pf+1]),reshape(plt_dif_qxx(:,:,sliceid),[pf+1,pf+1]))
end

colorbar
title('velocity x difference')

% difference in fluxes



% [rec.u.ux, rec.u.uy, rec.u.uz] = reconstruction.dim_3.surface(X(1:nr_surfc,1),rec.jacobian,p, rec.metric.invg, rec.basis.surfac);

% errors

% error_x = sqrt(sum((ext_qxx-rec.u.ux).^2 .* rec.metric.invg .* rec.wwwf))
% error_y = sqrt(sum((ext_qyy-rec.u.uy).^2 .* rec.metric.invg .* rec.wwwf))
% error_z = sqrt(sum((ext_qzz-rec.u.uz).^2 .* rec.metric.invg .* rec.wwwf))
