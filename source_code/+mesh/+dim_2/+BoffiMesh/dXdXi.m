function [new_dXdxii] = dXdXi(eleid, xii, eta, K, xbound)

% xbound = [0 pi];

meshDeformation = 0.45;
mySquare = mesh.dim_2.BoffiMesh.mimeticFEM2.IrregularSquare(K, meshDeformation);

mapping_factor = 0.5*(xbound(2)-xbound(1));

new_dXdxii = mapping_factor * mySquare.dXdXi(eleid, xii, eta);

end

