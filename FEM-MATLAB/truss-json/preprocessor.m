function  [K,f,d] = preprocessor(DataFile);
include_flags;

% input all variables 
FEMData = jsondecode(fileread(DataFile));

Title = FEMData.Title;
nsd  = FEMData.nsd;
ndof = FEMData.ndof;
nnp  = FEMData.nnp;
nel  = FEMData.nel;
nen  = FEMData.nen;

neq  = ndof*nnp;    % Number of equations

f 	= zeros(neq,1);   % Initialize force vector
d 	= zeros(neq,1);   % Initialize displacement matrix
K 	= zeros(neq);     % Initialize stiffness matrix

% Element properties
CArea = FEMData.CArea;
E     = FEMData.E;
 
% prescribed displacements
nd = FEMData.nd; 	% Number of prescribed displacement degrees-of-freedom
d  = FEMData.d;

% prescribed forces
fdof = FEMData.fdof;        % dofs in which forces are applied
f(fdof) = FEMData.force;    % forces applied in the dofs
 
% output plots
plot_truss 	= FEMData.plot_truss;
plot_node	= FEMData.plot_node;

x = FEMData.x;  % X coordinate  
y = FEMData.y;  % Y coordinate

% connectivity array
IEN = FEMData.IEN;     

leng = sqrt((x(IEN(:,2))-x(IEN(:,1))).^2 + (y(IEN(:,2))-y(IEN(:,1))).^2);

% plot truss
plottruss;

% Generate LM array 

for e = 1:nel
    for j = 1:nen
        for m = 1:ndof
            ind = (j-1)*ndof + m;
            LM(ind,e) = ndof*IEN(e, j) - ndof + m;
        end
    end
end
