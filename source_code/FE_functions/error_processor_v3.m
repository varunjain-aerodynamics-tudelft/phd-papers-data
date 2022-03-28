function [ error, l2err_phi ] = error_processor_v3(val_exact, val_calc, wfxy2, ggg)
%ERROR_PROCESSOR Summary of this function goes here
%   Detailed explanation goes here

% Created on 16 October, 2018
% Modified on 18 November, 2018

% this transpose the matrices along 3rd dimension 
% val_calc = permute(val_calc, [2 1 3]);

% check the transpose here (it should be on exact or calculated ?)
err_phi(:,:,:) = val_exact(:,:,:) - val_calc(:,:,:);

err_phi = (err_phi).^2 .* wfxy2 .* ggg;
l2err_phi = sum(sum(sum(err_phi)));
l2err_phi2 = sqrt(l2err_phi);

error.sqrt = l2err_phi2;
error.sqre = l2err_phi;

end

