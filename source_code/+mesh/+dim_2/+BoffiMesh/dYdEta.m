function [new_dYdeta] = dYdEta(eleid, xii, eta, K, ybound)

% ybound = [0 pi];

meshDeformation = 0.45;
mySquare = mesh.dim_2.BoffiMesh.mimeticFEM2.IrregularSquare(K, meshDeformation);

mapping_factor = 0.5*(ybound(2)-ybound(1));

new_dYdeta = mapping_factor * mySquare.dYdEta(eleid, xii, eta);

end

