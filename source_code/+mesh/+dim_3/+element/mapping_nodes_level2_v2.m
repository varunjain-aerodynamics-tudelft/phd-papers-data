function [xii,eta,zta] = mapping_nodes_level2_v2(s,t,u,el_bounds1,elx2,ely2,elz2,el_bounds2)

% xii1 = el_bounds1.x(elx1);
% xii2 = el_bounds1.x(elx1+1);
% 
% eta1 = el_bounds1.y(ely1);
% eta2 = el_bounds1.y(ely1+1);
% 
% zta1 = el_bounds1.z(elz1);
% zta2 = el_bounds1.z(elz1+1);

% new_el_bounds2.x = 0.5*(xii1 + xii2) + 0.5*(xii2 - xii1) * el_bounds2.x;
% new_el_bounds2.y = 0.5*(eta1 + eta2) + 0.5*(eta2 - eta1) * el_bounds2.y;
% new_el_bounds2.z = 0.5*(zta1 + zta2) + 0.5*(zta2 - zta1) * el_bounds2.z;

new_el_bounds2.x = 0.5*(el_bounds1.x(1:Kx_coarse) + el_bounds1.x(2:Kx_coarse+1)) + 0.5*(el_bounds1.x(2:Kx_coarse+1) - el_bounds1.x(1:Kx_coarse)) * el_bounds2.x;
new_el_bounds2.y = 0.5*(el_bounds1.y(1:Ky_coarse) + el_bounds1.y(2:Ky_coarse+1)) + 0.5*(el_bounds1.y(2:Ky_coarse+1) - el_bounds1.y(1:Ky_coarse)) * el_bounds2.y;
new_el_bounds2.z = 0.5*(el_bounds1.z(1:Kz_coarse) + el_bounds1.z(2:Kz_coarse+1)) + 0.5*(el_bounds1.z(2:Kz_coarse+1) - el_bounds1.z(1:Kz_coarse)) * el_bounds2.z;

new_xii1 = new_el_bounds2.x(elx2);
new_xii2 = new_el_bounds2.x(elx2+1);

new_eta1 = new_el_bounds2.y(ely2);
new_eta2 = new_el_bounds2.y(ely2+1);

new_zta1 = new_el_bounds2.z(elz2);
new_zta2 = new_el_bounds2.z(elz2+1);

xii = 0.5*(new_xii1 + new_xii2) + 0.5*(new_xii2 - new_xii1) * s;
eta = 0.5*(new_eta1 + new_eta2) + 0.5*(new_eta2 - new_eta1) * t;
zta = 0.5*(new_zta1 + new_zta2) + 0.5*(new_zta2 - new_zta1) * u;

end

