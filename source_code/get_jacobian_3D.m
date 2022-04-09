function [dtize_jacobian] = get_jacobian_3D(p,Kx,Ky,Kz,domain_jacobian)

ttl_nr_el = Kx*Ky*Kz;

el_bounds.x = linspace(-1,1,Kx+1);
el_bounds.y = linspace(-1,1,Ky+1);
el_bounds.z = linspace(-1,1,Kz+1);

element.nodes = @(xii,eta,zta,elx,ely,elz) mesh.dim_3.element.mapping_nodes(xii,eta,zta,elx,ely,elz,el_bounds);

%----------- create discretization nodes & basis

[x,w] = GLLnodes(p);
[h,e] = MimeticpolyVal(x,p,1);

h = eye(p+1); %%%%% this is only valid if you are using GLL Nodes - please check before using 

dtize.basis.hx = h;
dtize.basis.hy = h;
dtize.basis.hz = h;

dtize.basis.ex = e;
dtize.basis.ey = e;
dtize.basis.ez = e;

dtize.weights.wx = w;
dtize.weights.wy = w;
dtize.weights.wz = w;

dtize.nodes.sx = x;
dtize.nodes.sy = x;
dtize.nodes.sz = x;

%----------- create discretization jacobian and nodes

[dtize.nodes.eta, dtize.nodes.xii, dtize.nodes.zta] = meshgrid(dtize.nodes.sx,dtize.nodes.sy,dtize.nodes.sz);

dtize.nodes.xii = dtize.nodes.xii(:)';
dtize.nodes.eta = dtize.nodes.eta(:)';
dtize.nodes.zta = dtize.nodes.zta(:)';

for i = 1:Kx
    for j = 1:Ky
        for k = 1:Kz
            eleid = (k-1)*Kx*Ky + (j-1)*Kx + i;
            [dtize.nodes.xii2(:,:,eleid),dtize.nodes.eta2(:,:,eleid),dtize.nodes.zta2(:,:,eleid)] = element.nodes(dtize.nodes.xii, dtize.nodes.eta, dtize.nodes.zta,i,j,k);
        end
    end
end

% [dtize.nodes.x,dtize.nodes.y,dtize.nodes.z] = domain.mapping(dtize.nodes.xii2, dtize.nodes.eta2, dtize.nodes.zta2);

dtize.jacobian = domain_jacobian(dtize.nodes.xii2, dtize.nodes.eta2, dtize.nodes.zta2);

dXii_ds = 0.5*(el_bounds.x(2:end) - el_bounds.x(1:end-1));
dEta_dt = 0.5*(el_bounds.y(2:end) - el_bounds.y(1:end-1));
dZta_du = 0.5*(el_bounds.z(2:end) - el_bounds.z(1:end-1));

[d1, d2, d3] = meshgrid(dXii_ds,dEta_dt,dZta_du);

d1 = permute(d1,[2 1 3]);
d1 = d1(:)';

d2 = permute(d2,[2 1 3]);
d2 = d2(:)';

d3 = d3(:)';

dtize_jacobian.dXdxii = dtize.jacobian.dXdxii .* reshape(d1, [1 1 ttl_nr_el]);
dtize_jacobian.dYdxii = dtize.jacobian.dYdxii .* reshape(d1, [1 1 ttl_nr_el]);
dtize_jacobian.dZdxii = dtize.jacobian.dZdxii .* reshape(d1, [1 1 ttl_nr_el]);

dtize_jacobian.dXdeta = dtize.jacobian.dXdeta .* reshape(d2, [1 1 ttl_nr_el]);
dtize_jacobian.dYdeta = dtize.jacobian.dYdeta .* reshape(d2, [1 1 ttl_nr_el]);
dtize_jacobian.dZdeta = dtize.jacobian.dZdeta .* reshape(d2, [1 1 ttl_nr_el]);

dtize_jacobian.dXdzta = dtize.jacobian.dXdzta .* reshape(d3, [1 1 ttl_nr_el]);
dtize_jacobian.dYdzta = dtize.jacobian.dYdzta .* reshape(d3, [1 1 ttl_nr_el]);
dtize_jacobian.dZdzta = dtize.jacobian.dZdzta .* reshape(d3, [1 1 ttl_nr_el]);

end

