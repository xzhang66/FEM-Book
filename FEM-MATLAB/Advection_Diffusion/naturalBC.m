% compute and assemble nodal boundary force vector
function f = naturalBC(f);
Global_variables;

for i = 1:nnp
   if flags(i) == 1
      node = ID(i);
      f(node) = f(node) + Area(node)*n_bc(node);
   end
end