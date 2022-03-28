function [new_dXdeta] = dXdEta(eleid, xii, eta, K, ybound)

% ybound = [0 pi];

meshDeformation = 0.45;
mySquare = mesh.dim_2.BoffiMesh.mimeticFEM2.IrregularSquare(K, meshDeformation);

mapping_factor = 0.5*(ybound(2)-ybound(1));
% mapping_factor = 2/(ybound(2)-ybound(1));

new_dXdeta = mapping_factor * mySquare.dXdEta(eleid, xii, eta);

end

