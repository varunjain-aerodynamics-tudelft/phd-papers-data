function [dtize] = initialize_discretization_nodes_basis(p)

[x,w] = GLLnodes(p);
[h,e] = MimeticpolyVal(x,p,1);

dtize.basis.hx = h;
dtize.basis.hy = h;
dtize.basis.hz = h;

dtize.basis.ex = e;
dtize.basis.ey = e;
dtize.basis.ez = e;

dtize.weights.wx = w;
dtize.weights.wy = w;
dtize.weights.wz = w;

dtize.nodes.sx = x;
dtize.nodes.sy = x;
dtize.nodes.sz = x;

end

