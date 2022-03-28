function [ xi, eta ] = mapping_element_xii_eta( mapping, element_bounds_x, element_bounds_y, el_x, el_y, s, t )
%MAPPING_ELEMENT Summary of this function goes here
%   Detailed explanation goes here

% Created on : 8 March, 2019

% el = (el_x - 1)*K + el_y;

% element_bounds_x = linspace(-1, 1, K+1);
% element_bounds_y = linspace(-1, 1, K+1);

xi_i1  = element_bounds_x(el_x);
xi_i2  = element_bounds_x(el_x+1);
eta_j1 = element_bounds_y(el_y);
eta_j2 = element_bounds_y(el_y+1);

xi  = 0.5*(xi_i1 + xi_i2)   + 0.5*(xi_i2 - xi_i1).* s;
eta = 0.5*(eta_j1 + eta_j2) + 0.5*(eta_j2 - eta_j1).* t;

% [x, y] = mapping(xbound, ybound, c, xi ,eta);
% [x, y] = mapping(xi ,eta);

end

