function [domain] = quad_mapping(x1,x2,x3,x4,y1,y2,y3,y4)

lineqn = @(s,z1,z2) s*(z2-z1) + z1;
dlineqn_ds = @(s,z1,z2) ones(size(s)).*(z2-z1);

mapp.F_L_x = @(s,t) lineqn(t,x1,x4);
mapp.F_L_y = @(s,t) lineqn(t,y1,y4);

mapp.F_R_x = @(s,t) lineqn(t,x2,x3);
mapp.F_R_y = @(s,t) lineqn(t,y2,y3);

mapp.F_B_x = @(s,t) lineqn(s,x1,x2);
mapp.F_B_y = @(s,t) lineqn(s,y1,y2);

mapp.F_T_x = @(s,t) lineqn(s,x4,x3);
mapp.F_T_y = @(s,t) lineqn(s,y4,y3);

derr.dF_L_x = @(s,t) dlineqn_ds(t,x1,x4);
derr.dF_L_y = @(s,t) dlineqn_ds(t,y1,y4);

derr.dF_R_x = @(s,t) dlineqn_ds(t,x2,x3);
derr.dF_R_y = @(s,t) dlineqn_ds(t,y2,y3);

derr.dF_B_x = @(s,t) dlineqn_ds(s,x1,x2);
derr.dF_B_y = @(s,t) dlineqn_ds(s,y1,y2);

derr.dF_T_x = @(s,t) dlineqn_ds(s,x4,x3);
derr.dF_T_y = @(s,t) dlineqn_ds(s,y4,y3);

domain.mapping = @(xii,eta) mesh.transfinite_mesh.mapping(xii,eta,mapp);
domain.jacobian = @(xii,eta) mesh.transfinite_mesh.jacobian(xii,eta,mapp,derr);

end

