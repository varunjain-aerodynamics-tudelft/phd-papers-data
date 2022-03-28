function [dx_dxii] = dX_dxii(xii,eta,F_L_x,F_R_x,dFBx_ds,dFTx_ds)
%DX_DXII Summary of this function goes here
%   Detailed explanation goes here

% Created on 5 April, 2019

s = (xii+1)/2;
t = (eta+1)/2;

F00_x = F_L_x(0,0);
F0h_x = F_L_x(0,1);
Fh0_x = F_R_x(1,0);
F11_x = F_R_x(1,1);

dx_ds = - F_L_x(0,t)    + F_R_x(1,t)    + (1-t).*dFBx_ds(s,0)   + t.*dFTx_ds(s,1)...
        + (1-t)*F00_x   + t*F0h_x       - (1-t)*Fh0_x           - t*F11_x;

ds_dxii = 0.5;

dx_dxii = dx_ds * ds_dxii;

end

