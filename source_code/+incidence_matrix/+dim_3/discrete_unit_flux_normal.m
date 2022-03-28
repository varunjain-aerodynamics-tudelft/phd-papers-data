function [local_N] = discrete_unit_flux_normal(p,local_bndry_dof)

N.left = -1*ones(1,p^2);
N.rght = +1*ones(1,p^2);

N.bttm = -1*ones(1,p^2);
N.topp = +1*ones(1,p^2);

N.back = -1*ones(1,p^2);
N.frnt = +1*ones(1,p^2);

nr_local_lambda = 6*p^2;
nr_surfc = 3*p^2*(p+1);
nr_volum = p^3;
ttl_loc_dof = nr_surfc+nr_volum;

NN = [N.rght N.topp N.frnt N.left N.bttm N.back];

% local_bndry_dof = [globl_dof_rght globl_dof_topp globl_dof_frnt globl_dof_left globl_dof_bttm globl_dof_back];

local_N = sparse(1:nr_local_lambda,local_bndry_dof,NN,nr_local_lambda,ttl_loc_dof);

end

