function [ der ] = dx_deta( xbound, c, xi, eta )
%DX_DETA Summary of this function goes here
%   Detailed explanation goes here

% Created on : 22 March,2018

der = 0.5*(xbound(2)-xbound(1))*c*pi*sin(pi*xi).*cos(pi*eta);

end

