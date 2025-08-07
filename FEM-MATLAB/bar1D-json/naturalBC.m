% compute and assemble nodal boundary force vector
function f = naturalBC(f);
include_flags;

for i = 1:nnp
   if flags(i) == 1
      node = ID(i);
      f(node) = f(node) + CArea(i)*n_bc(i);
   end
end