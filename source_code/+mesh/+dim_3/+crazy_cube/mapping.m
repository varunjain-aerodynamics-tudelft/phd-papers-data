function [x,y,z] = mapping(xbound,ybound,zbound,c,xii,eta,zta)

x1 = xbound(1);
x2 = xbound(2);

y1 = ybound(1);
y2 = ybound(2);

z1 = zbound(1);
z2 = zbound(2);

x = (x1 + x2)/2 + (x2 - x1)/2 .* (xii + c .* sin(pi.*xii) .* sin(pi.*eta) .* sin(pi.*zta));
y = (y1 + y2)/2 + (y2 - y1)/2 .* (eta + c .* sin(pi.*xii) .* sin(pi.*eta) .* sin(pi.*zta));
z = (z1 + z2)/2 + (z2 - z1)/2 .* (zta + c .* sin(pi.*xii) .* sin(pi.*eta) .* sin(pi.*zta));

end

