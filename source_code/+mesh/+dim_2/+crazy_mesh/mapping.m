function [ x, y ] = mapping( x_bound, y_bound, c, xi, eta )
%MAPPING Summary of this function goes here
%   Detailed explanation goes here

% Created on : 22 March,2018

x1 = x_bound(1);
x2 = x_bound(2);
y1 = y_bound(1);
y2 = y_bound(2);

x = (x1 + x2)/2 + (x2 - x1)/2 .* (xi + c .* sin(pi .* xi) .* sin(pi .* eta));
y = (y1 + y2)/2 + (y2 - y1)/2 .* (eta + c .* sin(pi .* xi) .* sin(pi .* eta));

end

