clear all;
close all;

include_flags;

% Generate Input Data for Example 9.3 (64 element mesh)
FEMData.Title = "Example 9.3 (64 element mesh)";               
                   
% Mesh specifications
nsd = 2;          % number of spacial dimensions
nnp = 81;         % number of total nodes
nel = 64;         % number of total elements
nen = 4;          % number of nodes on each element
ndof = 2;         % degrees of freedom per node
neq  = nnp*ndof;  % number of equations

FEMData.nsd = nsd;
FEMData.ndof = ndof;
FEMData.nnp = nnp;
FEMData.nel = nel;
FEMData.nen = nen;

% material properties
E  = 30e6;     % Young's modulus 
ne = 0.3;      % Poisson ratio  

FEMData.E = E;
FEMData.ne = ne;

nd    = 10;            % number of essential boundary conditions
FEMData.nd = nd;

flags = zeros(neq,1);  % array to set B.C flags 
e_bc  = zeros(neq,1);  % essential B.C array
n_bc  = zeros(neq,1);  % natural B.C array

ngp   = 2;             % number of gauss points in each direction
FEMData.ngp = ngp;

% boundary conditions and point forces
nd = 18;     % number of essential boundary conditions
FEMData.nd = nd;

% essential B.C.
ind1 = 1:18:(73-1)*ndof+1;                 % all x dofs along the line y=0
ind2 = 2:18:(73-1)*ndof+2;                 % all y dofs along the line x=0 
flags(ind1) = 2;      e_bc(ind1)  = 0.0;
flags(ind2) = 2;      e_bc(ind2)  = 0.0;

FEMData.flags = flags;

% plots
plot_mesh      = 'yes';
plot_nod       = 'yes';
plot_disp      = 'yes';
compute_stress = 'yes';
plot_stress_xx = 'yes';
plot_mises     = 'yes';
fact           = 9.221e3;      % factor for scaled displacements plot

FEMData.plot_mesh = plot_mesh; 
FEMData.plot_nod = plot_nod; 
FEMData.plot_disp = plot_disp; 
FEMData.compute_stress = compute_stress; 
FEMData.plot_stress_xx = plot_stress_xx; 
FEMData.plot_mises = plot_mises; 
FEMData.fact = fact;       

% natural B.C  - given on edges 
n_bc    = [  73    74     75     76     77    78    79     80     % node 1
             74    75     76     77     78    79    80     81     % node 2 
             0      0      0      0      0     0     0      0     % traction at node 1 in x 
            -20   -20    -20    -20    -20   -20   -20    -20     % traction at node 1 in y
             0      0      0      0      0     0     0      0     % traction at node2  in x
            -20   -20    -20    -20    -20   -20   -20    -20];   % traction at node 2 in y
nbe   = 8;         

FEMData.n_bc = n_bc;
FEMData.nbe = nbe;

% mesh generation
mesh2d;

FEMData.x = x;
FEMData.y = y;
FEMData.IEN = IEN;

saveToJSON(FEMData, 'elasticity_64.json');