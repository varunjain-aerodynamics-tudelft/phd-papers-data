function [jacobian] = jacobian(xbound,ybound,zbound,xii,eta,zta)

x1 = xbound(1);
x2 = xbound(2);

y1 = ybound(1);
y2 = ybound(2);

z1 = zbound(1);
z2 = zbound(2);

cx = 0;
cy = 0;
cz = 0;

s = (x1 + x2)/2 + (x2 - x1)/2 .* (xii);
t = (y1 + y2)/2 + (y2 - y1)/2 .* (eta);
u = (z1 + z2)/2 + (z2 - z1)/2 .* (zta);

ds_dxii = (x2 - x1)/2;
dt_deta = (y2 - y1)/2;
du_dzta = (z2 - z1)/2;

jacobian.dXdxii = (1 - 3*pi*cx * sin(3*pi*s) .* cos(3*pi*t) .* cos(3*pi*u)) .* ds_dxii;
jacobian.dXdeta =    -(3*pi*cx * cos(3*pi*s) .* sin(3*pi*t) .* cos(3*pi*u)) .* dt_deta;
jacobian.dXdzta =    -(3*pi*cx * cos(3*pi*s) .* cos(3*pi*t) .* sin(3*pi*u)) .* du_dzta;

jacobian.dYdxii =    -(3*pi*cy * sin(3*pi*s) .* cos(3*pi*t) .* cos(3*pi*u)) .* ds_dxii;
jacobian.dYdeta = (1 - 3*pi*cy * cos(3*pi*s) .* sin(3*pi*t) .* cos(3*pi*u)) .* dt_deta;
jacobian.dYdzta =    -(3*pi*cy * cos(3*pi*s) .* cos(3*pi*t) .* sin(3*pi*u)) .* du_dzta;

jacobian.dZdxii =    -(3*pi*cz * sin(3*pi*s) .* cos(3*pi*t) .* cos(3*pi*u)) .* ds_dxii;
jacobian.dZdeta =    -(3*pi*cz * cos(3*pi*s) .* sin(3*pi*t) .* cos(3*pi*u)) .* dt_deta;
jacobian.dZdzta = (1 - 3*pi*cz * cos(3*pi*s) .* cos(3*pi*t) .* sin(3*pi*u)) .* du_dzta;

end

