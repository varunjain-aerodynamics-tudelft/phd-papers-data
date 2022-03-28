function [dtize_jacobian] = get_jacobian_2D(p,Kx,Ky,domain_jacobian)

ttl_nr_el = Kx*Ky;

el_bounds.x = linspace(-1,1,Kx+1);
el_bounds.y = linspace(-1,1,Ky+1);

element.nodes = @(xii,eta,elx,ely) mesh.dim_2.element.mapping_nodes(xii,eta,elx,ely,el_bounds);

%----------- create discretization nodes & basis

[x,w] = GLLnodes(p);
[h,e] = MimeticpolyVal(x,p,1);

% h = eye(p+1); %%%%% this is only valid if you are using GLL Nodes - please check before using 

dtize.basis.hx = h;
dtize.basis.hy = h;

dtize.basis.ex = e;
dtize.basis.ey = e;

dtize.weights.wx = w;
dtize.weights.wy = w;

dtize.nodes.sx = x;
dtize.nodes.sy = x;

%----------- create discretization jacobian and nodes

[dtize.nodes.eta, dtize.nodes.xii] = meshgrid(dtize.nodes.sx,dtize.nodes.sy);

dtize.nodes.xii = dtize.nodes.xii(:)';
dtize.nodes.eta = dtize.nodes.eta(:)';

for i = 1:Kx
    for j = 1:Ky
        eleid = (j-1)*Kx + i;
        [dtize.nodes.xii2(:,:,eleid),dtize.nodes.eta2(:,:,eleid)] = element.nodes(dtize.nodes.xii, dtize.nodes.eta,i,j);
    end
end

% [dtize.nodes.x,dtize.nodes.y,dtize.nodes.z] = domain.mapping(dtize.nodes.xii2, dtize.nodes.eta2, dtize.nodes.zta2);

dtize.jacobian = domain_jacobian(dtize.nodes.xii2, dtize.nodes.eta2);

dXii_ds = 0.5*(el_bounds.x(2:end) - el_bounds.x(1:end-1));
dEta_dt = 0.5*(el_bounds.y(2:end) - el_bounds.y(1:end-1));

[d1, d2] = meshgrid(dXii_ds,dEta_dt);

d1 = permute(d1,[2 1 3]);
d1 = d1(:)';

d2 = permute(d2,[2 1 3]);
d2 = d2(:)';

dtize_jacobian.dXdxii = dtize.jacobian.dXdxii .* reshape(d1, [1 1 ttl_nr_el]);
dtize_jacobian.dYdxii = dtize.jacobian.dYdxii .* reshape(d1, [1 1 ttl_nr_el]);

dtize_jacobian.dXdeta = dtize.jacobian.dXdeta .* reshape(d2, [1 1 ttl_nr_el]);
dtize_jacobian.dYdeta = dtize.jacobian.dYdeta .* reshape(d2, [1 1 ttl_nr_el]);

end