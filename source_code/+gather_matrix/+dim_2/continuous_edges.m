function [gcf] = continuous_edges(K,p)

% Created on : 3 May, 2020

gcf = zeros(2*p*(p+1), K^2);

B = 1:p*(p-1);
C = p*(p-1) + 1: 2*p*(p-1);

local_hor_edges = (p+1):p^2;
% local_ver_edges = p*(p-1) + 1: 2*p*(p-1);

for j = 1:p
    local_ver_edges((j-1)*(p-1) + (1:(p-1))) = p*(p+1) + (j-1)*(p+1) + (2:p);
end

for el = 1:K^2
    gcf(local_hor_edges,el) = B + 2*p*(p-1)*(el-1);
    gcf(local_ver_edges,el) = C + 2*p*(p-1)*(el-1);
end

% ggv = gather_matrix.GM_boundary_LM_nodes(K,p);
ggv = gather_matrix.dim_2.dual_surfaces_lambda(K,K,p);
ggv2 = ggv + K^2 * 2 * p * (p-1);

local_bttm = 1:p;
local_topp = p^2+1:p^2+p;
local_left = p*(p+1) + (1:p+1:p*(p+1));
local_rght = p*(p+1) + (p+1:p+1:p*(p+1)+p);

gcf(local_bttm,:) = ggv2(1:p,:);
gcf(local_topp,:) = ggv2(3*p+1:4*p,:);
gcf(local_left,:) = ggv2(2*p+1:3*p,:);
gcf(local_rght,:) = ggv2(p+1:2*p,:);

end

