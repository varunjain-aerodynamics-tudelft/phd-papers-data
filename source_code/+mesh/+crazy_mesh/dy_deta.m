function [ der ] = dy_deta( ybound, c, xi, eta )
%DY_DETA Summary of this function goes here
%   Detailed explanation goes here

% Created on : 22 March,2018

der = 0.5*(ybound(2)-ybound(1))*(1+ c*pi*sin(pi*xi).*cos(pi*eta));

end

