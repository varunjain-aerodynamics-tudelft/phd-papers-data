function [local_NK] = unit_normal_frnt_boundary(p, Kx, Ky, Kz, globl_dof_rght, globl_dof_left, globl_dof_topp, globl_dof_bttm, globl_dof_frnt, globl_dof_back)

%%%%% changed on 9 March, does not include P dof, only U's
ttl_nr_el = Kx * Ky * Kz;

nr_surfc = 3*p^2*(p+1);
nr_volum = p^3;
ttl_loc_dof = nr_surfc+nr_volum;

% gather_local_dof = gather_matrix.GM_total_local_dof_v3(ttl_nr_el, ttl_loc_dof);
% gather_local_dof = gather_local_dof';

[gather_qqq] = gather_matrix.dim_3.continuous_surfaces(Kx, Ky, Kz, p);

% NK.left = -1*ones(1,p^2*Ky*Kz);
% NK.rght = +1*ones(1,p^2*Ky*Kz);
% 
% NK.bttm = -1*ones(1,p^2*Kx*Kz);
% NK.topp = +1*ones(1,p^2*Kx*Kz);
% 
% NK.back = -1*ones(1,p^2*Kx*Ky);
NK.frnt = +1*ones(1,p^2*Kx*Ky);

% NNK = [NK.rght NK.left NK.topp NK.bttm NK.frnt NK.back];
% NNK = [NK.rght NK.topp NK.frnt NK.left NK.bttm NK.back];

% nr_boundary_lambda = 2*p^2*(Ky*Kz + Kx*Kz + Kx*Ky);
nr_boundary_lambda_frnt = p^2*Kx*Ky;

% this is for continuous elements
% ttl_loc_dof2 = p^2 * (Ky*Kz)* (Kx*p+1) + p^2 * (Kx*Kz)* (Ky*p+1) + p^2 * (Kx*Ky)* (Kz*p+1);

% ttl_loc_dof2 = Kx*p*Ky*p*(Kz*p+1) + Ky*p*Kz*p*(Kx*p+1) + Kx*p*Kz*p*(Ky*p+1) + Kx*Ky*Kz*p^3;
ttl_loc_dof3 = Kx*p*Ky*p*(Kz*p+1) + Ky*p*Kz*p*(Kx*p+1) + Kx*p*Kz*p*(Ky*p+1);

% % right face 
% for j = 1:Ky
%     for k = 1:Kz
%         local_eleid = (k-1)*Ky + j;
%         globl_eleid = (k-1)*Kx*Ky + (j-1)*Kx + Kx;
%         boundary_face_rght(:,local_eleid) = gather_qqq(globl_dof_rght,globl_eleid); %% gather local dof will change as required for continuous elements 
%     end
% end
% 
% boundary_face_rght = boundary_face_rght(:)';
% 
% % left face
% for j = 1:Ky
%     for k = 1:Kz
%         local_eleid = (k-1)*Ky + j;
%         globl_eleid = (k-1)*Kx*Ky + (j-1)*Kx + 1;
%         boundary_face_left(:,local_eleid) = gather_qqq(globl_dof_left,globl_eleid);
%     end
% end
% 
% boundary_face_left = boundary_face_left(:)';
% 
% % top face
% for i = 1:Kx
%     for k = 1:Kz
%         local_eleid = (k-1)*Kx + i;
%         globl_eleid = (k-1)*Kx*Ky + (Ky-1)*Kx + i;
%         boundary_face_topp(:,local_eleid) = gather_qqq(globl_dof_topp,globl_eleid);
%     end
% end
% 
% boundary_face_topp = boundary_face_topp(:)';
% 
% % bottom face
% for i = 1:Kx
%     for k = 1:Kz
%         local_eleid = (k-1)*Kx + i;
%         globl_eleid = (k-1)*Kx*Ky + (1-1)*Kx + i;
%         boundary_face_bttm(:,local_eleid) = gather_qqq(globl_dof_bttm,globl_eleid);
%     end
% end
% 
% boundary_face_bttm = boundary_face_bttm(:)';

% front face 
for i = 1:Kx
    for j = 1:Ky
        local_eleid = (j-1)*Kx + i;
        globl_eleid = (Kz-1)*Kx*Ky + (j-1)*Kx + i;
        boundary_face_frnt(:,local_eleid) = gather_qqq(globl_dof_frnt,globl_eleid);
    end
end

boundary_face_frnt = boundary_face_frnt(:)';

% % back face 
% for i = 1:Kx
%     for j = 1:Ky
%         local_eleid = (j-1)*Kx + i;
%         globl_eleid = (1-1)*Kx*Ky + (j-1)*Kx + i;
%         boundary_face_back(:,local_eleid) = gather_qqq(globl_dof_back,globl_eleid);
%     end
% end
% 
% boundary_face_back = boundary_face_back(:)';

% local_bndry_dof_NK = [boundary_face_rght boundary_face_left boundary_face_topp boundary_face_bttm boundary_face_frnt boundary_face_back];
% local_bndry_dof_NK = [boundary_face_rght boundary_face_topp boundary_face_frnt boundary_face_left boundary_face_bttm boundary_face_back];

local_NK = sparse(1:nr_boundary_lambda_frnt,boundary_face_frnt,NK.frnt,nr_boundary_lambda_frnt,ttl_loc_dof3);

end

