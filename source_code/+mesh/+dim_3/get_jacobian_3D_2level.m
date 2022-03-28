function [dtize] = get_jacobian_3D_2level(dtize,Kx_coarse,Ky_coarse,Kz_coarse,Kx_fine,Ky_fine,Kz_fine,domain)

ttl_nr_sub_el = Kx_fine*Ky_fine*Kz_fine;

ttl_nr_mac_el = Kx_coarse*Ky_coarse*Kz_coarse;

ttl_nr_el = ttl_nr_sub_el * ttl_nr_mac_el;

el_bounds1.x = linspace(-1,1,Kx_coarse+1);
el_bounds1.y = linspace(-1,1,Ky_coarse+1);
el_bounds1.z = linspace(-1,1,Kz_coarse+1);

el_bounds2.x = linspace(-1,1,Kx_fine+1);
el_bounds2.y = linspace(-1,1,Ky_fine+1);
el_bounds2.z = linspace(-1,1,Kz_fine+1);

element.nodes2 = @(xii,eta,zta,elx1,ely1,elz1,elx2,ely2,elz2) mesh.dim_3.element.mapping_nodes_level2(xii,eta,zta,elx1,ely1,elz1,el_bounds1,elx2,ely2,elz2,el_bounds2);

% here eta is first to fix the anomoly of meshgrid function 
[dtize.nodes.eta, dtize.nodes.xii, dtize.nodes.zta] = meshgrid(dtize.nodes.sx,dtize.nodes.sy,dtize.nodes.sz);

dtize.nodes.xii = dtize.nodes.xii(:)';
dtize.nodes.eta = dtize.nodes.eta(:)';
dtize.nodes.zta = dtize.nodes.zta(:)';

dtize.www = kron(kron(dtize.weights.wz,dtize.weights.wy),dtize.weights.wx);

for i1 = 1:Kx_coarse
    for j1 = 1:Ky_coarse
        for k1 = 1:Kz_coarse
            eleid_coarse = (k1-1)*Kx_coarse*Ky_coarse + (j1-1)*Kx_coarse + i1;
            for i2 = 1:Kx_fine
                for j2 = 1:Ky_fine
                    for k2 = 1:Kz_fine
                        eleid_fine = (eleid_coarse-1)*ttl_nr_sub_el + (k2-1)*Kx_fine*Ky_fine + (j2-1)*Kx_fine + i2;
                        [dtize.nodes.xii2(:,:,eleid_fine),dtize.nodes.eta2(:,:,eleid_fine),dtize.nodes.zta2(:,:,eleid_fine)] = element.nodes2(dtize.nodes.xii, dtize.nodes.eta, dtize.nodes.zta,i1,j1,k1,i2,j2,k2);
                    end
                end
            end
        end
    end
end

[dtize.nodes.x,dtize.nodes.y,dtize.nodes.z] = domain.mapping(dtize.nodes.xii2,dtize.nodes.eta2,dtize.nodes.zta2);

dtize.jacobian = domain.jacobian(dtize.nodes.xii2,dtize.nodes.eta2,dtize.nodes.zta2);

dXii_ds1 = 0.5*(el_bounds1.x(2:end) - el_bounds1.x(1:end-1));
dEta_dt1 = 0.5*(el_bounds1.y(2:end) - el_bounds1.y(1:end-1));
dZta_du1 = 0.5*(el_bounds1.z(2:end) - el_bounds1.z(1:end-1));

[d1_1, d1_2, d1_3] = meshgrid(dXii_ds1,dEta_dt1,dZta_du1);

% d1_1 = permute(d1_1,[2 1 3]);
d1_1 = d1_1(:)';

% d1_2 = permute(d1_2,[2 1 3]);
d1_2 = d1_2(:)';

d1_3 = d1_3(:)';

dXii_ds2 = 0.5*(el_bounds2.x(2:end) - el_bounds2.x(1:end-1));
dEta_dt2 = 0.5*(el_bounds2.y(2:end) - el_bounds2.y(1:end-1));
dZta_du2 = 0.5*(el_bounds2.z(2:end) - el_bounds2.z(1:end-1));

[d2_1, d2_2, d2_3] = meshgrid(dXii_ds2,dEta_dt2,dZta_du2);

% d2_1 = permute(d2_1,[2 1 3]);
d2_1 = d2_1(:)';

% d2_2 = permute(d2_2,[2 1 3]);
d2_2 = d2_2(:)';

d2_3 = d2_3(:)';

d1 = kron(d1_1,d2_1);
d2 = kron(d1_2,d2_2);
d3 = kron(d1_3,d2_3);

dtize.jacobian.dXdxii = dtize.jacobian.dXdxii .* reshape(d1, [1 1 ttl_nr_el]);
dtize.jacobian.dYdxii = dtize.jacobian.dYdxii .* reshape(d1, [1 1 ttl_nr_el]);
dtize.jacobian.dZdxii = dtize.jacobian.dZdxii .* reshape(d1, [1 1 ttl_nr_el]);

dtize.jacobian.dXdeta = dtize.jacobian.dXdeta .* reshape(d2, [1 1 ttl_nr_el]);
dtize.jacobian.dYdeta = dtize.jacobian.dYdeta .* reshape(d2, [1 1 ttl_nr_el]);
dtize.jacobian.dZdeta = dtize.jacobian.dZdeta .* reshape(d2, [1 1 ttl_nr_el]);

dtize.jacobian.dXdzta = dtize.jacobian.dXdzta .* reshape(d3, [1 1 ttl_nr_el]);
dtize.jacobian.dYdzta = dtize.jacobian.dYdzta .* reshape(d3, [1 1 ttl_nr_el]);
dtize.jacobian.dZdzta = dtize.jacobian.dZdzta .* reshape(d3, [1 1 ttl_nr_el]);

end

