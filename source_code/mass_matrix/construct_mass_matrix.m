function [mass_matrix] = construct_mass_matrix(test_func,trial_func,weights,metric)

trial_func = trial_func .* weights .* metric;

mass_matrix = test_func * trial_func';

end

