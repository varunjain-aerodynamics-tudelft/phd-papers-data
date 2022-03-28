function [jacobian] = jacobian(xbound,ybound,c,xii,eta)

x1 = xbound(1);
x2 = xbound(2);

y1 = ybound(1);
y2 = ybound(2);

jacobian.dXdxii = 0.5*(x2-x1)*(1+ c*pi*cos(pi*xii).*sin(pi*eta));
jacobian.dXdeta = 0.5*(x2-x1)*c*pi*sin(pi*xii).*cos(pi*eta);

jacobian.dYdxii = 0.5*(y2-y1)*c*pi*cos(pi*xii).*sin(pi*eta);
jacobian.dYdeta = 0.5*(y2-y1)*(1+ c*pi*sin(pi*xii).*cos(pi*eta));

end

