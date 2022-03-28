function [cochain] = reduction(ff_int, dual_basis, weights, metric)

ff_int = ff_int .* weights .* metric;

cochain = dual_basis * ff_int';

end

