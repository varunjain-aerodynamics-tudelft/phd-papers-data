function [x,y,z] = mapping(xbound,ybound,xii,eta)

x1 = xbound(1);
x2 = xbound(2);

y1 = ybound(1);
y2 = ybound(2);

cx = 0.03;
cy = -0.04;

s = (x1 + x2)/2 + (x2 - x1)/2 .* (xii);
t = (y1 + y2)/2 + (y2 - y1)/2 .* (eta);

x = s + cx .* cos(3*pi.*s) .* cos(3*pi.*t);
y = t + cy .* cos(3*pi.*s) .* cos(3*pi.*t);

end