function [g] = poisson_metric(jac,invg)

g.g11 = (jac.dXdxii .* jac.dXdxii + jac.dYdxii .* jac.dYdxii + jac.dZdxii .* jac.dZdxii) .* invg;
g.g12 = (jac.dXdxii .* jac.dXdeta + jac.dYdxii .* jac.dYdeta + jac.dZdxii .* jac.dZdeta) .* invg;
g.g13 = (jac.dXdxii .* jac.dXdzta + jac.dYdxii .* jac.dYdzta + jac.dZdxii .* jac.dZdzta) .* invg;

g.g21 = (jac.dXdeta .* jac.dXdxii + jac.dYdeta .* jac.dYdxii + jac.dZdeta .* jac.dZdxii) .* invg;
g.g22 = (jac.dXdeta .* jac.dXdeta + jac.dYdeta .* jac.dYdeta + jac.dZdeta .* jac.dZdeta) .* invg;
g.g23 = (jac.dXdeta .* jac.dXdzta + jac.dYdeta .* jac.dYdzta + jac.dZdeta .* jac.dZdzta) .* invg;

g.g31 = (jac.dXdzta .* jac.dXdxii + jac.dYdzta .* jac.dYdxii + jac.dZdzta .* jac.dZdxii) .* invg;
g.g32 = (jac.dXdzta .* jac.dXdeta + jac.dYdzta .* jac.dYdeta + jac.dZdzta .* jac.dZdeta) .* invg;
g.g33 = (jac.dXdzta .* jac.dXdzta + jac.dYdzta .* jac.dYdzta + jac.dZdzta .* jac.dZdzta) .* invg;

end

