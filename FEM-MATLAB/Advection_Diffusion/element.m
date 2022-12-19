function [ ke, fe ] = element( e, a )
%ELEMENT Element matrix computations

Global_variables;

IENe = IEN(:,e);          % extract local connectivity information
xe = x(IENe);             % extract element x coordinates
lene = xe(nen) - xe(1);   % the length of the element
J = lene/2;               % compute Jacobian
[w , gp] = gauss(ngp);    % extract Gauss points and weights

ke = zeros(ndof*nen,ndof*nen);      % initialize element stiffness matrix
fe = zeros(ndof*nen,1);             % initialize element nodal force vector

for i = 1:ngp
   xt = 0.5*(xe(1)+xe(nen))+J*gp(i);  % Compute Gauss points in physical coordinates

   N = Nmatrix1D(xt,xe);    % shape functions matrix
   B = Bmatrix1D(xt,xe);    % derivative of shape functions matrix

   Ae = N*Area(IENe);       % cross-sectional area at element gauss points
   kce = N*k(IENe);         % Diffusion coefficient at element gauss points
   be = N*body(IENe);       % body forces at element gauss points
   ke = ke + w(i)*J*(B'*Ae*(kce*(1+PN*a))*B)... % Diffusion matrix is modified by alpha
       + w(i)*J*(N'*Ae*(-2*kce*PN/lene)*B);     % Advection matrix
   fe = fe + w(i)*N'*be;         % compute element nodal body heat source vector
end

end

