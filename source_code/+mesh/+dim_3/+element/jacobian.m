function [jacobian_new] = jacobian(s,t,u,elx,ely,elz,el_bounds,domain_jacobian)

xii1 = el_bounds.x(elx);
xii2 = el_bounds.x(elx+1);

eta1 = el_bounds.y(ely);
eta2 = el_bounds.y(ely+1);

zta1 = el_bounds.z(elz);
zta2 = el_bounds.z(elz+1);

xii = 0.5*(xii1 + xii2) + 0.5*(xii2 - xii1) * s;
eta = 0.5*(eta1 + eta2) + 0.5*(eta2 - eta1) * t;
zta = 0.5*(zta1 + zta2) + 0.5*(zta2 - zta1) * u;

dxii_ds = 0.5*(el_bounds.x(elx+1) - el_bounds.x(elx));
deta_dt = 0.5*(el_bounds.y(ely+1) - el_bounds.y(ely));
dzta_du = 0.5*(el_bounds.z(elz+1) - el_bounds.z(elz));

jacobian_old = domain_jacobian(xii,eta,zta);

jacobian_new.dXdxii = jacobian_old.dXdxii * dxii_ds;
jacobian_new.dXdeta = jacobian_old.dXdeta * deta_dt;
jacobian_new.dXdzta = jacobian_old.dXdzta * dzta_du;

jacobian_new.dYdxii = jacobian_old.dYdxii * dxii_ds;
jacobian_new.dYdeta = jacobian_old.dYdeta * deta_dt;
jacobian_new.dYdzta = jacobian_old.dYdzta * dzta_du;

jacobian_new.dZdxii = jacobian_old.dZdxii * dxii_ds;
jacobian_new.dZdeta = jacobian_old.dZdeta * deta_dt;
jacobian_new.dZdzta = jacobian_old.dZdzta * dzta_du;

end

