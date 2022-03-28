function [gather_qqq] = continuous_surfaces(Kx, Ky, Kz, p)

nr_surfc = 3*p^2*(p+1);
ttl_nr_el = Kx * Ky * Kz;
local_bndry_face = mesh.dim_3.element_boundary_faces_map(p);

gather_qqq = zeros(nr_surfc,ttl_nr_el);

% yz faces / right faces 
count = 1;
for i = 2:p
    for j = 1:p
        for k = 1:p
            local_face_yz(count) = (k-1)*p*(p+1) + (j-1)*(p+1) + i;
            globl_face_yz(count) = (k-1)*p*(p-1) + (j-1)*(p-1) + (i-1);
            count = count+1;
        end
    end
end

% xz faces / bottom faces
count = 1;
for i = 1:p
    for j = 2:p
        for k = 1:p
            local_face_xz(count) = p^2*(p+1) + (k-1)*p*(p+1) + (j-1)*p + i;
            globl_face_xz(count) = p^2*(p-1) + (k-1)*p*(p-1) + (j-1 -1)*p + i;
            count = count+1;
        end
    end
end

%  xy faces / front faces 
count = 1;
for i = 1:p
    for j = 1:p
        for k = 2:p
            local_face_xy(count) = 2*p^2*(p+1) + (k-1)*p^2 + (j-1)*p + i;
            globl_face_xy(count) = 2*p^2*(p-1) + (k-1-1)*p^2 + (j-1)*p + i;
            count = count+1;
        end
    end
end

% for elx = 1:Kx
%     for ely = 1:Ky
%         for elz = 1:Kz
%             eln = (elz-1)*Kx*Ky + (ely-1)*Kx + elx;
local_face_rght = local_bndry_face.rght;
local_face_topp = local_bndry_face.topp;
local_face_frnt = local_bndry_face.frnt;
local_face_left = local_bndry_face.left;
local_face_bttm = local_bndry_face.bttm;
local_face_back = local_bndry_face.back;

bndry_faces = [local_face_rght local_face_topp local_face_frnt local_face_left local_face_bttm local_face_back];
%         end
%     end
% end

gather_lambda = gather_matrix.dim_3.lambda_pressure(Kx,Ky,Kz,p);

% for elx = 1:Kx
%     for ely = 1:Ky
%         for elz = 1:Kz
%             eln = (elz-1)*Kx*Ky + (ely-1)*Kx + elx;
%             if p~=1
%                 gather_qqq(local_face_yz, eln) = (eln-1)*3*p^2*(p-1) + globl_face_yz;
%                 gather_qqq(local_face_xz, eln) = (eln-1)*3*p^2*(p-1) + globl_face_xz;
%                 gather_qqq(local_face_xy, eln) = (eln-1)*3*p^2*(p-1) + globl_face_xy;
%             end
%             gather_qqq(bndry_faces,eln) = ttl_nr_el*3*p^2*(p-1) + gather_lambda(:,eln);
%         end
%     end
% end

for eln = 1:ttl_nr_el
    if p~=1
        gather_qqq(local_face_yz, eln) = (eln-1)*3*p^2*(p-1) + globl_face_yz;
        gather_qqq(local_face_xz, eln) = (eln-1)*3*p^2*(p-1) + globl_face_xz;
        gather_qqq(local_face_xy, eln) = (eln-1)*3*p^2*(p-1) + globl_face_xy;
    end
    gather_qqq(bndry_faces,eln) = ttl_nr_el*3*p^2*(p-1) + gather_lambda(:,eln);
end

end

