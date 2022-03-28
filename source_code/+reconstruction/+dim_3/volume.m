function [reconstruction] = volume(dof,basis,pf,ggg)

basis = basis ./ ggg;

reconstruction = basis' * dof;

% reconstruction = reshape(reconstruction, [pf+1 pf+1 pf+1]);

end

