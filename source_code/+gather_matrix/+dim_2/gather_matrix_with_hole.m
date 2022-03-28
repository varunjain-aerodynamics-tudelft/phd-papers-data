function [gather_matrix] = gather_matrix_with_hole(gather_matrix,Kx,Ky,p)

A = 1:p;

idx_bttm_lambda = A;
idx_rght_lambda = p + A;
idx_left_lambda = 2*p + A;
idx_topp_lambda = 3*p + A;

ELEID = reshape(1:Kx*Ky, [Kx Ky]);

% overlap edges

for j = 1:Ky
    eleid_rght = ELEID(Kx,j);
    gather_matrix(idx_rght_lambda, eleid_rght) = A + Kx*(Ky-1)*p + Ky*(Kx-1)*p + (j-1)*p;
    
    eleid_left = ELEID(1,j);
    gather_matrix(idx_left_lambda, eleid_left) = A + Kx*(Ky-1)*p + Ky*(Kx-1)*p + (j-1)*p;
end

% bndry edges

for i = 1:Kx
    eleid_bot = ELEID(i,1);
    gather_matrix(idx_bttm_lambda, eleid_bot) = A + Kx*(Ky-1)*p + Ky*(Kx-1)*p + Ky*p + (i-1)*p;
    
    eleid_topp = ELEID(i,Ky);
    gather_matrix(idx_topp_lambda, eleid_topp) = A + Kx*(Ky-1)*p + Ky*(Kx-1)*p + Ky*p + Kx*p + (i-1)*p;
end

end

