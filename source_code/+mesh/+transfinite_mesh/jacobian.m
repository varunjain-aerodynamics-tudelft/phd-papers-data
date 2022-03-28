function [jacobian_der] = jacobian(xii,eta,mapp,derr)

s = (xii+1)/2;
t = (eta+1)/2;

ds_dxii = 0.5;
dt_deta = 0.5;

F00_x = mapp.F_L_x(0,0);
F0h_x = mapp.F_L_x(0,1);
Fh0_x = mapp.F_R_x(1,0);
F11_x = mapp.F_R_x(1,1);

F00_y = mapp.F_B_y(0,0);
F0h_y = mapp.F_T_y(0,1);
Fh0_y = mapp.F_B_y(1,0);
F11_y = mapp.F_T_y(1,1);

%% ---------- dx_dxii

dx_ds = - mapp.F_L_x(0,t)    + mapp.F_R_x(1,t)    + (1-t).*derr.dF_B_x(s,0)   + t.*derr.dF_T_x(s,1)...
        + (1-t)*F00_x   + t*F0h_x       - (1-t)*Fh0_x           - t*F11_x;

%% ---------- dy_deta

dy_dt = (1-s).*derr.dF_L_y(0,t) + s.*derr.dF_R_y(1,t)   - mapp.F_B_y(s,0)     + mapp.F_T_y(s,1)...
        + (1-s)*F00_y       - (1-s)*F0h_y       + s*Fh0_y        - s*F11_y;

%% ---------- dx_deta

dx_dt = (1-s).*derr.dF_L_x(0,t) + s.*derr.dF_R_x(1,t)   - mapp.F_B_x(s,0)     + mapp.F_T_x(s,1)...
        + (1-s)*F00_x       - (1-s)*F0h_x       + s*Fh0_x        - s*F11_x;

%% ---------- dy_dxii

dy_ds = - mapp.F_L_y(0,t)    + mapp.F_R_y(1,t)    + (1-t).*derr.dF_B_y(s,0)   + t.*derr.dF_T_y(s,1)...
        + (1-t)*F00_y   + t*F0h_y       - (1-t)*Fh0_y           - t*F11_y;

%% ---------- map from [0,1] to [-1,1]

jacobian_der.dXdxii = dx_ds * ds_dxii;
jacobian_der.dYdeta = dy_dt * dt_deta;
jacobian_der.dXdeta = dx_dt * dt_deta;
jacobian_der.dYdxii = dy_ds * ds_dxii;

end

