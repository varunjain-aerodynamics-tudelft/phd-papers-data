function [dtize_metric] = get_metric_3D_darcy(dtize_jacobian,perm)

dtize_metric_ggg = metric.dim_3.determinant_jacobian(dtize_jacobian);
dtize_metric_invg = 1 ./dtize_metric_ggg;
dtize_metric = metric.dim_3.darcy_metric(dtize_jacobian, dtize_metric_invg, perm);

end

