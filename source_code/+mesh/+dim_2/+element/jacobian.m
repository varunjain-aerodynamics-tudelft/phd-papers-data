function [jacobian_new] = jacobian(s,t,elx,ely,el_bounds,domain_jacobian)

xii1 = el_bounds.x(elx);
xii2 = el_bounds.x(elx+1);

eta1 = el_bounds.y(ely);
eta2 = el_bounds.y(ely+1);

xii = 0.5*(xii1 + xii2) + 0.5*(xii2 - xii1) * s;
eta = 0.5*(eta1 + eta2) + 0.5*(eta2 - eta1) * t;

dxii_ds = 0.5*(el_bounds.x(elx+1) - el_bounds.x(elx));
deta_dt = 0.5*(el_bounds.y(ely+1) - el_bounds.y(ely));

jacobian_old = domain_jacobian(xii,eta);

jacobian_new.dXdxii = jacobian_old.dXdxii * dxii_ds;
jacobian_new.dXdeta = jacobian_old.dXdeta * deta_dt;

jacobian_new.dYdxii = jacobian_old.dYdxii * dxii_ds;
jacobian_new.dYdeta = jacobian_old.dYdeta * deta_dt;

end

