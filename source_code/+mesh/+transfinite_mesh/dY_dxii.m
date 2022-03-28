function [dy_dxii] = dY_dxii(xii,eta,F_L_y,F_R_y,dFBy_ds,dFTy_ds)
%DY_DXII Summary of this function goes here
%   Detailed explanation goes here

% Created on 5 April, 2019

s = (xii+1)/2;
t = (eta+1)/2;

F00_y = F_L_y(0,0);
F0h_y = F_L_y(0,1);
Fh0_y = F_R_y(1,0);
F11_y = F_R_y(1,1);

dy_ds = - F_L_y(0,t)    + F_R_y(1,t)    + (1-t).*dFBy_ds(s,0)   + t.*dFTy_ds(s,1)...
        + (1-t)*F00_y   + t*F0h_y       - (1-t)*Fh0_y           - t*F11_y;

ds_dxii = 0.5;

dy_dxii = dy_ds * ds_dxii;
    
end

