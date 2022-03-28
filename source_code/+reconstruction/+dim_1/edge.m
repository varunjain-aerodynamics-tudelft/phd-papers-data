function [lambda_rec] = edge(primal_lambda_dof,ef)

p = size(primal_lambda_dof,1);
int_pts = size(ef,2);

lambda_rec = zeros(1,int_pts);

for i = 1:p
    for rec_pt = 1:int_pts
        lambda_rec(rec_pt) = lambda_rec(rec_pt) + primal_lambda_dof(i)*ef(i,rec_pt);
    end
end

lambda_rec = - lambda_rec;

end

