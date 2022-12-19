% Input File for Pelect Number 0.1

nsd = 1;     % number of space dimensions  
ndof =1;     % number of degrees-of-freedom per node 
nnp = 21;     % number of nodal points 
nel = 20;     % number of elements 
nen = 2;     % number of element nodes 

neq = ndof*nnp;         % number of equations 

f = zeros(neq,1);      % initialize nodal force vector 
d = zeros(neq,1);      % initialize nodal displacement vector 
K = zeros(neq);        % initialize stiffness matrix 

flags = zeros(neq,1);  % initialize flag vector 
e_bc = zeros(neq,1);   % initialize vector of essential boundary conditions 
n_bc = zeros(neq,1);   % initialize vector of natural boundary conditions 

% element characteristics (given at element nodes) 
k      =  [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]';      % Diffusion coefficient 
body   =  [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]';      % body forces 
Area   =  [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]';      % cross-sectional area 
PN     =  100;                      % Pelect Number
alpha(1)  =  0;                     % Galerkin method
alpha(2)  =  1;                     % One-Side FD
alpha(3)  =  coth(abs(PN))-1./abs(PN); % Petrov-Galerkin with optimal value for alpha

% gauss integration 
ngp    = 2;                          % number of gauss points 

% essential B.C.'s (displacements or temperatures) 
flags(1) = 2;                 % flags to mark nodes located for the essential boundary 
e_bc(1) = 0;                  % value of essential B.C 
flags(2) = 2;                % flags to mark nodes located for the essential boundary 
e_bc(2) = 1;                 % value of essential B.C 
nd      = 2;                  % number of nodes on the essential boundary 

% natural B.C.'s (stresses or fluxes) 
% there is no natural B.C. in this problem

% point forces or point heat source applied at any point 
% there is no point heat source in this problem
np = 0;

% mesh generation 
mesh_Advection_Diffusion;