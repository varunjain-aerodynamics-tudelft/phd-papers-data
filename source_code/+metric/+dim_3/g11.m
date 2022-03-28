function [g] = g11(jac, invg)

dXdxii = jac.dXdxii;
dYdxii = jac.dYdxii;
dZdxii = jac.dZdxii;

dXdeta = jac.dXdeta;
dYdeta = jac.dYdeta;
dZdeta = jac.dZdeta;

dXdzta = jac.dXdzta;
dYdzta = jac.dYdzta;
dZdzta = jac.dZdzta;

alpha_1 = (((dXdxii.*dZdzta - dXdzta.*dZdxii) .* (dXdxii.*dYdeta - dXdeta.*dYdxii))...
          - ((dXdxii.*dZdeta - dXdeta.*dZdxii) .* (dXdxii.*dYdzta - dXdzta.*dYdxii))) .* invg .* invg;
      
alpha_2 = (((dXdeta.*dZdzta - dXdzta.*dZdeta) .* (dXdxii.*dYdeta - dXdeta.*dYdxii))...
          - ((dXdxii.*dZdeta - dXdeta.*dZdxii) .* (dXdeta.*dYdzta) - dXdzta.*dYdeta)) .* invg .* invg;

alpha_3 = (((dXdeta.*dZdzta - dXdzta.*dZdeta) .* (dXdxii.*dYdzta - dXdzta.*dYdxii))...
          - ((dXdxii.*dZdzta - dXdzta.*dZdxii) .* (dXdeta.*dYdzta - dXdzta.*dYdeta))) .* invg .* invg;

beta_1 = (((dXdxii.*dYdeta - dXdeta.*dYdxii) .* (dYdxii.*dZdzta - dYdzta.*dZdxii))...
          - ((dXdxii.*dYdzta - dXdzta.*dYdxii) .* (dYdxii.*dZdeta - dYdeta.*dZdxii))) .* invg .* invg;
      
beta_2 = (((dXdxii.*dYdeta - dXdeta.*dYdxii) .* (dYdeta.*dZdzta - dYdzta.*dZdeta))...
          - ((dXdeta.*dYdzta - dXdzta.*dYdeta) .* (dYdxii.*dZdeta - dYdeta.*dZdxii))) .* invg .* invg;
      
beta_3 = (((dXdxii.*dYdzta - dXdzta.*dYdxii) .* (dYdeta.*dZdzta - dYdzta.*dZdeta))...
          - ((dXdeta.*dYdzta - dXdzta.*dYdeta) .* (dYdxii.*dZdzta - dYdzta.*dZdxii))) .* invg .* invg;

gamma_1 = ((dYdxii.*dZdzta - dYdzta.*dZdxii) .* (dXdxii .* dZdeta - dXdeta.*dZdxii)...
          - (dYdxii.*dZdeta - dYdeta.*dZdxii) .* (dXdxii.*dZdzta - dXdzta.*dZdxii)) .* invg .* invg;

gamma_2 = ((dYdeta.*dZdzta - dYdzta.*dZdeta) .* (dXdxii.*dZdeta - dXdeta.*dZdxii)...
          - (dYdxii.*dZdeta - dYdeta.*dZdxii) .* (dXdeta.*dZdzta - dXdzta.*dZdeta)) .* invg .* invg;

gamma_3 = ((dYdeta.*dZdzta - dYdzta.*dZdeta) .* (dXdxii.*dZdzta - dXdzta.*dZdxii)...
          - (dYdxii.*dZdzta - dYdzta.*dZdxii) .* (dXdeta.*dZdzta - dXdzta.*dZdeta)) .* invg .* invg;

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

