function [ gf ] = volumes( K, p )
%GM_2_FORM Summary of this function goes here
%   Detailed explanation goes here

% Created on : 19 Mar, 2018

gf = zeros(p^2, K^2);

el_total_dof = p^2;
num_el_dof_q = p^2;

q = 1:p^2;

for ele = 1:K^2
    gf(1:num_el_dof_q,ele) = q + el_total_dof*(ele-1);
end

end

