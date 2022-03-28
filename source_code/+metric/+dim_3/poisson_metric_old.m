function [g] = poisson_metric_old(jac,invg)

dXdxii = jac.dXdxii;
dYdxii = jac.dYdxii;
dZdxii = jac.dZdxii;

dXdeta = jac.dXdeta;
dYdeta = jac.dYdeta;
dZdeta = jac.dZdeta;

dXdzta = jac.dXdzta;
dYdzta = jac.dYdzta;
dZdzta = jac.dZdzta;

dXiidx = (jac.dYdeta.*jac.dZdzta - jac.dYdzta.*jac.dZdeta) .* invg;
dXiidy = -(jac.dYdxii.*jac.dZdzta - jac.dYdzta.*jac.dZdxii) .* invg;
dXiidz = (jac.dYdxii.*jac.dZdeta - jac.dYdeta.*jac.dZdxii) .* invg;

dEtadx = -(jac.dXdeta.*jac.dZdzta - jac.dXdzta.*jac.dZdeta) .* invg;
dEtady = (jac.dXdxii.*jac.dZdzta - jac.dXdzta.*jac.dZdxii) .* invg;
dEtadz = -(jac.dXdxii.*jac.dZdeta - jac.dXdeta.*jac.dZdxii) .* invg;

dZtadx = (jac.dXdeta.*jac.dYdzta - jac.dXdzta.*jac.dYdeta) .* invg;
dZtady = -(jac.dXdxii.*jac.dYdzta - jac.dXdzta.*jac.dYdxii) .* invg;
dZtadz = (jac.dXdxii.*jac.dYdeta - jac.dXdeta.*jac.dYdxii) .* invg;

alpha_1 = dEtady.*dZtadz - dEtadz.*dZtady;
alpha_2 = dEtadz.*dZtadx - dEtadx.*dZtadz;
alpha_3 = dEtadx.*dZtady - dEtady.*dZtadx;

beta_1 = dZtady.*dXiidz - dZtadz.*dXiidy;
beta_2 = dZtadz.*dXiidx - dZtadx.*dXiidz;
beta_3 = dZtadx.*dXiidy - dZtady.*dXiidx;

gamma_1 = dXiidy.*dEtadz - dXiidz.*dEtady;
gamma_2 = dXiidz.*dEtadx - dXiidx.*dEtadz;
gamma_3 = dXiidx.*dEtady - dXiidy.*dEtadx;

g.g11 = (alpha_1 .* dXdxii + alpha_2 .* dYdxii + alpha_3 .* dZdxii);
g.g12 = (beta_1  .* dXdxii + beta_2  .* dYdxii + beta_3  .* dZdxii);
g.g13 = (gamma_1 .* dXdxii + gamma_2 .* dYdxii + gamma_3 .* dZdxii);

g.g21 = (alpha_1 .* dXdeta + alpha_2 .* dYdeta + alpha_3 .* dZdeta);
g.g22 = (beta_1  .* dXdeta + beta_2  .* dYdeta + beta_3  .* dZdeta);
g.g23 = (gamma_1 .* dXdeta + gamma_2 .* dYdeta + gamma_3 .* dZdeta);

g.g31 = (alpha_1 .* dXdzta + alpha_2 .* dYdzta + alpha_3 .* dZdzta);
g.g32 = (beta_1  .* dXdzta + beta_2  .* dYdzta + beta_3  .* dZdzta);
g.g33 = (gamma_1 .* dXdzta + gamma_2 .* dYdzta + gamma_3 .* dZdzta);

end

