function [dy_deta] = dY_deta(xii,eta,F_B_y,F_T_y,dFLy_dt,dFRy_dt)
%DY_DETA Summary of this function goes here
%   Detailed explanation goes here

% Created on 5 April, 2019

s = (xii+1)/2;
t = (eta+1)/2;

F00_y = F_B_y(0,0);
F0h_y = F_T_y(0,1);
Fh0_y = F_B_y(1,0);
F11_y = F_T_y(1,1);

dy_dt = (1-s).*dFLy_dt(0,t) + s.*dFRy_dt(1,t)   - F_B_y(s,0)     + F_T_y(s,1)...
        + (1-s)*F00_y       - (1-s)*F0h_y       + s*Fh0_y        - s*F11_y;

dt_deta = 0.5;

dy_deta = dy_dt * dt_deta;

end

