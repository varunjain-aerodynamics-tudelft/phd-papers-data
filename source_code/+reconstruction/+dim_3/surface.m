function[ux,uy,uz] = surface(qh,jac,p,invg, basis)

qxh = qh(1:p^2*(p+1));
qyh = qh(p^2*(p+1)+1:2*p^2*(p+1));
qzh = qh(2*p^2*(p+1)+1:3*p^2*(p+1));

dXiidx = (jac.dYdeta.*jac.dZdzta - jac.dYdzta.*jac.dZdeta) .* invg;
dXiidy = -(jac.dYdxii.*jac.dZdzta - jac.dYdzta.*jac.dZdxii) .* invg;
dXiidz = (jac.dYdxii.*jac.dZdeta - jac.dYdeta.*jac.dZdxii) .* invg;

dEtadx = -(jac.dXdeta.*jac.dZdzta - jac.dXdzta.*jac.dZdeta) .* invg;
dEtady = (jac.dXdxii.*jac.dZdzta - jac.dXdzta.*jac.dZdxii) .* invg;
dEtadz = -(jac.dXdxii.*jac.dZdeta - jac.dXdeta.*jac.dZdxii) .* invg;

dZtadx = (jac.dXdeta.*jac.dYdzta - jac.dXdzta.*jac.dYdeta) .* invg;
dZtady = -(jac.dXdxii.*jac.dYdzta - jac.dXdzta.*jac.dYdxii) .* invg;
dZtadz = (jac.dXdxii.*jac.dYdeta - jac.dXdeta.*jac.dYdxii) .* invg;

% ------ velocity - x
basis_qx.yz = basis.yz .* (dEtady.*dZtadz - dEtadz.*dZtady);
basis_qx.zx = basis.zx .* (dZtady.*dXiidz - dZtadz.*dXiidy);
basis_qx.xy = basis.xy .* (dXiidy.*dEtadz - dXiidz.*dEtady);

ux1 = qxh' * basis_qx.yz;
ux2 = qyh' * basis_qx.zx;
ux3 = qzh' * basis_qx.xy;

ux = ux1 + ux2 + ux3;

% ------ velocity - y
basis_qy.yz = basis.yz .* (dEtadz.*dZtadx - dEtadx.*dZtadz);
basis_qy.zx = basis.zx .* (dZtadz.*dXiidx - dZtadx.*dXiidz);
basis_qy.xy = basis.xy .* (dXiidz.*dEtadx - dXiidx.*dEtadz);

uy1 = qxh' * basis_qy.yz;
uy2 = qyh' * basis_qy.zx;
uy3 = qzh' * basis_qy.xy;

uy = uy1 + uy2 + uy3;

% ------ velocity - z
basis_qz.yz = basis.yz .* (dEtadx.*dZtady - dEtady.*dZtadx);
basis_qz.zx = basis.zx .* (dZtadx.*dXiidy - dZtady.*dXiidx);
basis_qz.xy = basis.xy .* (dXiidx.*dEtady - dXiidy.*dEtadx);

uz1 = qxh' * basis_qz.yz;
uz2 = qyh' * basis_qz.zx;
uz3 = qzh' * basis_qz.xy;

uz = uz1 + uz2 + uz3;

end

