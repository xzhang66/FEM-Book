% File to include global variables 
global Title   % Title describing the problem to be solved
global nsd     % number of space dimensions
global ndof    % degrees of freedom per node 
global nnp     % number of nodal nodes
global nel     % number of elements
global nen     % number of element nodes
global neq     % number of equations
global ngp     % number of gauss points in each direction
global nd      % number of essential boundary conditions (x and y)
global e_bc    % essential B.C array
global np      % apply nodal forces ? (1: yes, 0: no)
global P       % point forces applied at nodes 
global b       % body forces defined at nodes
global D       % elasticity matrix
global LM      % location matrix
global ID      % identification array
global IEN     % connectivity array
global flags   % boundary condition flags
global n_bc    % natural B.C array
global x       % X coordinates 
global y       % Y coordinates
global nbe     % number of edges on the boundary

global counter         % counter of nodes for stress plots 
global nodestress      % stresses at nodes for the stress plots [sxx syy sxy] 
global compute_stress  % computer stresses ?

global plot_mesh       % plot mesh ?
global plot_disp       % plot displaced mesh ?
global plot_nod        % plot node number ?
global plot_stress_xx  % plot stress xx ?
global plot_mises      % plot mises stress ?
global fact            % factor for scaled displacements plot
