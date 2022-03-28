function [x,y,z] = mapping(xbound,ybound,zbound,xii,eta,zta)

x1 = xbound(1);
x2 = xbound(2);

y1 = ybound(1);
y2 = ybound(2);

z1 = zbound(1);
z2 = zbound(2);

cx = 0.03;
cy = -0.04;
cz = 0.05;

s = (x1 + x2)/2 + (x2 - x1)/2 .* (xii);
t = (y1 + y2)/2 + (y2 - y1)/2 .* (eta);
u = (z1 + z2)/2 + (z2 - z1)/2 .* (zta);

x = s + cx .* cos(3*pi.*s) .* cos(3*pi.*t) .* cos(3*pi.*u);
y = t + cy .* cos(3*pi.*s) .* cos(3*pi.*t) .* cos(3*pi.*u);
z = u + cz .* cos(3*pi.*s) .* cos(3*pi.*t) .* cos(3*pi.*u);

end

