function [dx_deta] = dX_deta(xii,eta,F_B_x,F_T_x,dFLx_dt,dFRx_dt)
%DX_DETA Summary of this function goes here
%   Detailed explanation goes here

% Created on 5 April, 2019

s = (xii+1)/2;
t = (eta+1)/2;

F00_x = F_B_x(0,0);
F0h_x = F_T_x(0,1);
Fh0_x = F_B_x(1,0);
F11_x = F_T_x(1,1);

dx_dt = (1-s).*dFLx_dt(0,t) + s.*dFRx_dt(1,t)   - F_B_x(s,0)     + F_T_x(s,1)...
        + (1-s)*F00_x       - (1-s)*F0h_x       + s*Fh0_x        - s*F11_x;

dt_deta = 0.5;

dx_deta = dx_dt * dt_deta;

end

