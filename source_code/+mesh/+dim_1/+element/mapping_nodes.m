function [xii] = mapping_nodes(s,elx,el_bounds)

xii1 = el_bounds.x(elx);
xii2 = el_bounds.x(elx+1);

xii = 0.5*(xii1 + xii2) + 0.5*(xii2 - xii1) * s;

end

