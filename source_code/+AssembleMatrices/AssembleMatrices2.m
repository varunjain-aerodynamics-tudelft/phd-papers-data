function [ C2 ] = AssembleMatrices2( GMrow, GMcol, data )

%ASSEMBLEMATRICES Summary of this function goes here
%   Detailed explanation goes here

% Created on  : 8th Feb, 2018
% Modified on : 10th Feb, 2018 

[loc_row, loc_col, num_el] = size(data);

cc = zeros(loc_row, loc_col, num_el);
rr = zeros(loc_row, loc_col, num_el);

tot_row = max(max(GMrow));
tot_col = max(max(GMcol));

for i = 1:num_el
    r = GMrow(:,i);
    c = GMcol(:,i);

    [cc(:,:,i), rr(:,:,i)] = meshgrid(c,r);
end

cc   = cc(:);
rr   = rr(:);
data = data(:);

C2 = sparse(rr, cc, data, tot_row, tot_col);

end