function [jacobian] = jacobian(xbound,ybound,zbound,xii,eta,zta)

x1 = xbound(1);
x2 = xbound(2);

y1 = ybound(1);
y2 = ybound(2);

z1 = zbound(1);
z2 = zbound(2);

cx = 0.03;
cy = -0.04;
cz = 0.05;

jacobian.dXdxii = (x2-x1)/2 * (1 - 3*pi*cx * sin(3*pi*xii) .* cos(3*pi*eta) .* cos(3*pi*zta));
jacobian.dXdeta = -(x2-x1)/2 * (3*pi*cx * cos(3*pi*xii) .* sin(3*pi*eta) .* cos(3*pi*zta));
jacobian.dXdzta = -(x2-x1)/2 * (3*pi*cx * cos(3*pi*xii) .* cos(3*pi*eta) .* sin(3*pi*zta));

jacobian.dYdxii = -(y2-y1)/2 * (3*pi*cy * sin(3*pi*xii) .* cos(3*pi*eta) .* cos(3*pi*zta));
jacobian.dYdeta = (y2-y1)/2 * (1 - 3*pi*cy * cos(3*pi*xii) .* sin(3*pi*eta) .* cos(3*pi*zta));
jacobian.dYdzta = -(y2-y1)/2 * (3*pi*cy * cos(3*pi*xii) .* cos(3*pi*eta) .* sin(3*pi*zta));

jacobian.dZdxii = -(z2-z1)/2 * (3*pi*cz * sin(3*pi*xii) .* cos(3*pi*eta) .* cos(3*pi*zta));
jacobian.dZdeta = -(z2-z1)/2 * (3*pi*cz * cos(3*pi*xii) .* sin(3*pi*eta) .* cos(3*pi*zta));
jacobian.dZdzta = (z2-z1)/2 * (1 - 3*pi*cz * cos(3*pi*xii) .* cos(3*pi*eta) .* sin(3*pi*zta));

end

