function [xii,eta] = mapping_nodes(s,t,elx,ely,el_bounds)

xii1 = el_bounds.x(elx);
xii2 = el_bounds.x(elx+1);

eta1 = el_bounds.y(ely);
eta2 = el_bounds.y(ely+1);

xii = 0.5*(xii1 + xii2) + 0.5*(xii2 - xii1) * s;
eta = 0.5*(eta1 + eta2) + 0.5*(eta2 - eta1) * t;

end

