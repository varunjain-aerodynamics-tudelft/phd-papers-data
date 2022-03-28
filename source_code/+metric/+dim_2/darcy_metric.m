function [metric] = darcy_metric(jac,perm)

dx_dxi  = jac.dXdxii;
dx_deta = jac.dXdeta;
dy_dxi  = jac.dYdxii;
dy_deta = jac.dYdeta;

ggg = jac.dXdxii .* jac.dYdeta - jac.dXdeta .* jac.dYdxii;

metric.invg = 1./ggg;

K11 = perm.kxx;
K12 = perm.kxy;
K22 = perm.kyy;

detK = K11 .* K22 - K12.^2;

k11i =  K22./detK;
k12i = -K12./detK;
k22i =  K11./detK;

metric.g11 = (k11i.*dx_deta.^2 + 2.*k12i.*dx_deta.*dy_deta + k22i.*dy_deta.^2)./ggg;
metric.g12 = (k11i.*dx_dxi.*dx_deta + k12i.*(dy_dxi.*dx_deta + dx_dxi.*dy_deta) + k22i.*dy_dxi.*dy_deta)./ggg;
metric.g22 = (k11i.*dx_dxi.^2 + 2.* k12i.*dy_dxi.*dx_dxi + k22i.*dy_dxi.^2)./ggg;

end

