function [local_bndry_face] = element_boundary_faces_map(p)

local_bndry_face.left = zeros(1,p^2);
local_bndry_face.rght = zeros(1,p^2);
local_bndry_face.bttm = zeros(1,p^2);
local_bndry_face.topp = zeros(1,p^2);
local_bndry_face.back = zeros(1,p^2);
local_bndry_face.frnt = zeros(1,p^2);

% left face
for j = 1:p
    for k = 1:p
        local_dof = (k-1)*p + j;
        local_bndry_face.left(local_dof) = (k-1)*p*(p+1) + (j-1)*(p+1) + 1;
    end
end

% right face
for j = 1:p
    for k = 1:p
        local_dof = (k-1)*p + j;
        local_bndry_face.rght(local_dof) = (k-1)*p*(p+1) + (j-1)*(p+1) + p + 1;
    end
end

% bottom face

for i = 1:p
    for k = 1:p
        local_dof = (k-1)*p + i;
        local_bndry_face.bttm(local_dof) = p^2*(p+1) + (k-1)*p*(p+1) + i;
    end
end

% top face

for i = 1:p
    for k = 1:p
        local_dof = (k-1)*p + i;
        local_bndry_face.topp(local_dof) = p^2*(p+1) + (k-1)*p*(p+1) + i + p^2;
    end
end

% back face

for i = 1:p
    for j = 1:p
        local_dof = (j-1)*p+i;
        local_bndry_face.back(local_dof) = 2*p^2*(p+1) + local_dof;
    end
end

% front face

for i = 1:p
    for j = 1:p
        local_dof = (j-1)*p+i;
        local_bndry_face.frnt(local_dof) = 2*p^2*(p+1) + local_dof + p^3;
    end
end

end

