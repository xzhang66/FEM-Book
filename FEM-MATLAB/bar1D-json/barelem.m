% generate element stiffness matrix and element nodal body force vector
function [ke, fe] = barelem(e);
include_flags;

IENe = IEN(:,e);          % extract local connectivity information
xe = x(IENe);             % extract element x coordinates
J = (xe(nen) - xe(1))/2;  % compute Jacobian
[w , gp] = gauss(ngp);    % extract Gauss points and weights

ke = zeros(nen,nen);      % initialize element stiffness matrix
fe = zeros(nen,1);        % initialize element nodal force vector

for i = 1:ngp
   xt = 0.5*(xe(1)+xe(nen))+J*gp(i);  % Compute Gauss points in physical coordinates

   N = Nmatrix1D(xt,xe);    % shape functions matrix
   B = Bmatrix1D(xt,xe);    % derivative of shape functions matrix

   Ae = N*CArea(IENe);      % cross-sectional area at element gauss points
   Ee = N*E(IENe);          % Young's modulus at element gauss points
   be = N*body(IENe);       % body forces at element gauss points
   ke = ke + w(i)*(B'*Ae*Ee*B);  % compute element stiffness matrix
   fe = fe + w(i)*N'*be;         % compute element nodal body force vector
end

ke = J*ke;
fe = J*fe;

% check for point forces in this element
for i=1:np        % loop over all point forces
   Pi = P(i);     % extract point force
   xpi = xp(i);   % extract the location of point force within an element
   if xe(1)<=xpi & xpi<xe(nen)
      fe = fe + Pi*[Nmatrix1D(xpi,xe)]';   % add to the nodal force vector
   end
end