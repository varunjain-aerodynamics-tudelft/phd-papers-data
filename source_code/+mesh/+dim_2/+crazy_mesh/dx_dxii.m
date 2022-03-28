function [ der ] = dx_dxii( xbound, c, xi, eta )
%DX_DXI Summary of this function goes here
%   Detailed explanation goes here

% Created on : 22 March,2018

der = 0.5*(xbound(2)-xbound(1))*(1+ c*pi*cos(pi*xi).*sin(pi*eta));

end

