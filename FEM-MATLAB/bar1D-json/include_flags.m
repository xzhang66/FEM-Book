% Include global variables
global Title   % Title describing the problem to be solved
global nsd  % number of space dimensions
global ndof % number of degrees-of-freedom per node
global nnp  % number of nodal points
global nel  % number of elements
global nen  % number of element nodes
global neq  % number of equations
global nd   % number of nodes on the essential boundary
global CArea  % nodal values of cross-sectional area
global E    % nodal values Young's modulus

global flags    % flags to mark nodes located on a boundary
                %   = 2: located on the essential boundary
                %   = 1: located on the natural boundary
global ID 
global IEN  % Element connectivity array
global LM   % Location matrix
global body % nodal values body forces
global x    % x coordinate
global y    % y is used only for the bar plot

global np   % number of point forces
global xp   % array of coordinates where point forces are applied
global P    % array of point forcess

global ngp  % number of gauss points

global n_bc % value of natural B.C
global e_bc % value of essential B.C

global plot_bar 
global plot_nod 
global nplot    % number of points in the element to plot displacements and stresses


