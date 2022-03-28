function [ der ] = dy_dxii( ybound, c, xi, eta )
%DY_DXI Summary of this function goes here
%   Detailed explanation goes here

% Created on : 22 March,2018

der = 0.5*(ybound(2)-ybound(1))*c*pi*cos(pi*xi).*sin(pi*eta);

end

