function [jacobian] = jacobian(xbound,ybound,zbound,c,xii,eta,zta)

x1 = xbound(1);
x2 = xbound(2);

y1 = ybound(1);
y2 = ybound(2);

z1 = zbound(1);
z2 = zbound(2);

jacobian.dXdxii = (x2-x1)/2 * (1 + pi*c * cos(pi*xii) .* sin(pi*eta) .* sin(pi*zta));
jacobian.dXdeta = (x2-x1)/2 * (pi*c * sin(pi*xii) .* cos(pi*eta) .* sin(pi*zta));
jacobian.dXdzta = (x2-x1)/2 * (pi*c * sin(pi*xii) .* sin(pi*eta) .* cos(pi*zta));

jacobian.dYdxii = (y2-y1)/2 * (pi*c * cos(pi*xii) .* sin(pi*eta) .* sin(pi*zta));
jacobian.dYdeta = (y2-y1)/2 * (1 + pi*c * sin(pi*xii) .* cos(pi*eta) .* sin(pi*zta));
jacobian.dYdzta = (y2-y1)/2 * (pi*c * sin(pi*xii) .* sin(pi*eta) .* cos(pi*zta));

jacobian.dZdxii = (z2-z1)/2 * (pi*c * cos(pi*xii) .* sin(pi*eta) .* sin(pi*zta));
jacobian.dZdeta = (z2-z1)/2 * (pi*c * sin(pi*xii) .* cos(pi*eta) .* sin(pi*zta));
jacobian.dZdzta = (z2-z1)/2 * (1 + pi*c * sin(pi*xii) .* sin(pi*eta) .* cos(pi*zta));

end

