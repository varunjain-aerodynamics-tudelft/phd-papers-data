function [xii2,eta2,zta2] = new_ref_nodes(range_i,range_j,range_z,xii,eta,zta,Ky,element_nodes)

% eleid = @(i,j,k) (k-1)*Kx*Ky + (j-1)*Kx + i;

left_eleid = @(i,j,k) (k-1)*Ky + j;

% eleid = @(i,j,k) (k-1)*Ky + j;
% 
% eleid = @(i,j,k) (k-1)*Ky + j;

for i = range_i
    for j = range_j
        for k = range_z
            eleid = left_eleid(i,j,k);
            [xii2(:,:,eleid),eta2(:,:,eleid),zta2(:,:,eleid)] = element_nodes(xii,eta,zta,i,j,k);
        end
    end
end

end

