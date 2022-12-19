function [K,f,d] = preprocessor(DataFile); 
include_flags;

% input all variables 
FEMData = jsondecode(fileread(DataFile));

Title = FEMData.Title;
nsd  = FEMData.nsd;
ndof = FEMData.ndof;
nnp  = FEMData.nnp;
nel  = FEMData.nel;
nen  = FEMData.nen;

neq  = ndof*nnp;      % Number of equations

f 	= zeros(neq,1);   % Initialize force vector
d 	= zeros(neq,1);   % Initialize displacement matrix
K 	= zeros(neq);     % Initialize stiffness matrix

% material properties
E  = FEMData.E;
ne = FEMData.ne;

D  = E/(1-ne^2) * [1    ne     0           
                   ne    1     0   
                   0     0     (1-ne)/2]; 
 
counter    = zeros(nnp,1);  % counter of nodes for stress plots 
nodestress = zeros(nnp,3);  % stresses at nodes for the stress plots [sxx syy sxy] 
 
flags = zeros(neq,1);  % array to set B.C flags  
e_bc  = zeros(neq,1);  % essential B.C array 
 
P     = zeros(neq,1);          % point forces applied at nodes 
b     = zeros(nen*ndof,nel);   % body forces defined at nodes 

np    = FEMData.np;
if np > 0
    P = FEMData.P;
end

ngp = FEMData.ngp;
nd  = FEMData.nd;

flags = FEMData.flags;

plot_mesh      = FEMData.plot_mesh; 
plot_nod       = FEMData.plot_nod; 
plot_disp      = FEMData.plot_disp; 
compute_stress = FEMData.compute_stress; 
plot_stress_xx = FEMData.plot_stress_xx; 
plot_mises     = FEMData.plot_mises; 
fact           = FEMData.fact;       

n_bc = FEMData.n_bc;
nbe  = FEMData.nbe;
x    = transpose(FEMData.x);
y    = transpose(FEMData.y);
IEN  = FEMData.IEN;

plotmesh; 

% generate ID and LM arrays
d = setup_ID_LM(d);
