function [x] = chol_sol(A,b)

R = chol(A);
x = R\(R'\b);

end

