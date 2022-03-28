function [x,y] = mapping(xii,eta,mapp)

s = (xii+1)/2;
t = (eta+1)/2;

F00_x = mapp.F_L_x(0,0);
F0h_x = mapp.F_L_x(0,1);
Fh0_x = mapp.F_R_x(1,0);
F11_x = mapp.F_R_x(1,1);

F00_y = mapp.F_L_y(0,0);
F0h_y = mapp.F_L_y(0,1);
Fh0_y = mapp.F_R_y(1,0);
F11_y = mapp.F_R_y(1,1);

x = (1-s).*mapp.F_L_x(0,t)     + s.*mapp.F_R_x(1,t)   + (1-t).*mapp.F_B_x(s,0)     + t.*mapp.F_T_x(s,1)...
    - (1-s).*(1-t)*F00_x  - (1-s).*t*F0h_x  - (1-t).*s*Fh0_x        - t.*s*F11_x;

y = (1-s).*mapp.F_L_y(0,t)     + s.*mapp.F_R_y(1,t)   + (1-t).*mapp.F_B_y(s,0)     + t.*mapp.F_T_y(s,1)...
    - (1-s).*(1-t)*F00_y  - (1-s).*t*F0h_y  - (1-t).*s*Fh0_y        - t.*s*F11_y;

end