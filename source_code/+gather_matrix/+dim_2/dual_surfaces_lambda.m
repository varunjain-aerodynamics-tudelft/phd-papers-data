function [gather_matrix] = dual_surfaces_lambda(Kx,Ky,p)

num_el = Kx * Ky;
num_local_boundary = 4*p;

gather_matrix = zeros(num_local_boundary, num_el);

A = 1:p;

idx_bttm_lambda = A;
idx_rght_lambda = p + A;
idx_left_lambda = 2*p + A;
idx_topp_lambda = 3*p + A;

ELEID = reshape(1:Kx*Ky, [Kx Ky]);

% bottom edges
count = 0;
for i = 1:Kx
    for j = 2:Ky
        eln = ELEID(i,j);
        gather_matrix(idx_bttm_lambda, eln) = p*count + A;
        count = count+1;
    end
end

% left edges
count = 0;
for i = 2:Kx
    for j = 1:Ky
        eln = ELEID(i,j);
        gather_matrix(idx_left_lambda, eln) = Kx*(Ky-1)*p + count*p + A;
        count = count +1;
    end
end

% right edges
for i = 1:Kx-1
    for j = 1:Ky
        eln = ELEID(i,j);
        eln_rght = ELEID(i+1,j);
        gather_matrix(idx_rght_lambda, eln) = gather_matrix(idx_left_lambda, eln_rght);
    end
end

% top edges
for i = 1:Kx
    for j = 1:Ky-1
        eln = ELEID(i,j);
        eln_topp = ELEID(i,j+1);
        gather_matrix(idx_topp_lambda, eln) = gather_matrix(idx_bttm_lambda, eln_topp);
    end
end

% bndry edges

for i = 1:Kx
    eleid_bot = ELEID(i,1);
    gather_matrix(idx_bttm_lambda, eleid_bot) = A + Kx*(Ky-1)*p + Ky*(Kx-1)*p + (i-1)*p;
    
    eleid_topp = ELEID(i,Ky);
    gather_matrix(idx_topp_lambda, eleid_topp) = A + Kx*(Ky-1)*p + Ky*(Kx-1)*p + Kx*p + 2*Ky*p + (i-1)*p;
end

% bdry edges
for j = 1:Ky
    eleid_rght = ELEID(Kx,j);
    gather_matrix(idx_rght_lambda, eleid_rght) = A + Kx*(Ky-1)*p + Ky*(Kx-1)*p + Kx*p + (j-1)*p;
    
    eleid_left = ELEID(1,j);
    gather_matrix(idx_left_lambda, eleid_left) = A + Kx*(Ky-1)*p + Ky*(Kx-1)*p + Kx*p + Ky*p + (j-1)*p;
end

end

