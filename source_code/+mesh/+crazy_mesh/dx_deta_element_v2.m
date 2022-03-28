function [ der ] = dx_deta_element_v2(dx_deta, element_bounds_x, element_bounds_y, elx, ely, s, t, xi2, eta2)
%DX_DETA_ELEMENT Summary of this function goes here
%   Detailed explanation goes here

% Created on : 8 March, 2019

% xi_i1  = element_bounds_x(elx);
% xi_i2  = element_bounds_x(elx+1);
% eta_j1 = element_bounds_y(ely);
% eta_j2 = element_bounds_y(ely+1);

% xi  = 0.5*(xi_i1 + xi_i2)   + 0.5*(xi_i2 - xi_i1).* s;
% eta = 0.5*(eta_j1 + eta_j2) + 0.5*(eta_j2 - eta_j1).* t;

deta_dt = 0.5*(element_bounds_y(ely+1) - element_bounds_y(ely));

der = dx_deta(xi2 ,eta2) * deta_dt;

end

