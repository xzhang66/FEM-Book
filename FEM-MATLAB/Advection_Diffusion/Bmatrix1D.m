% derivative of the shape functions computed in the physical coordinate - xt
function B = Bmatrix1D(xt,xe)
Global_variables;

if nen == 2       % derivative of linear shape functions (constant)
   B = 1/(xe(1)-xe(2))*[-1 1];
elseif nen == 3   % derivative of quadratic shape functions
   B(1)=(2*xt-xe(2)-xe(3))/((xe(1)-xe(2))*(xe(1)-xe(3)));
   B(2)=(2*xt-xe(1)-xe(3))/((xe(2)-xe(1))*(xe(2)-xe(3)));
   B(3)=(2*xt-xe(1)-xe(2))/((xe(3)-xe(1))*(xe(3)-xe(2)));
end