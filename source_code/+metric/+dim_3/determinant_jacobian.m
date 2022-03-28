function [det] = determinant_jacobian(jac)

det =   jac.dXdxii.*(jac.dYdeta.*jac.dZdzta - jac.dYdzta.*jac.dZdeta)...
      - jac.dXdeta.*(jac.dYdxii.*jac.dZdzta - jac.dYdzta.*jac.dZdxii)...
      + jac.dXdzta.*(jac.dYdxii.*jac.dZdeta - jac.dYdeta.*jac.dZdxii);

end

