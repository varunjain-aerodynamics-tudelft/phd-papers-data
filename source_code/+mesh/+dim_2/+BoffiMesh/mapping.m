function [x2,y2] = mapping(eleid, xii, eta, K, xbound, ybound)

% xbound = [0 pi];
% ybound = [0 pi];

meshDeformation = 0.45;
mySquare = mesh.dim_2.BoffiMesh.mimeticFEM2.IrregularSquare(K, meshDeformation);

[x,y] = mySquare.mapping(eleid, xii, eta);

x2 = 0.5*(xbound(2)+xbound(1)) + 0.5*x*(xbound(2)-xbound(1));
y2 = 0.5*(ybound(2)+ybound(1)) + 0.5*y*(ybound(2)-ybound(1));

end

