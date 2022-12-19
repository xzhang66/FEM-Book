% preprocessing: reads input data and sets up mesh information
%   Originally by Haim Waisman, Rensselaer
%   Modified by Xiong Zhang to read input data from json file
%
function [K,f,d] = preprocessor(DataFile);
include_flags;

% input all variables 
FEMData = jsondecode(fileread(DataFile));

Title = FEMData.Title;
nsd  = FEMData.nsd;     % number of space dimensions
ndof = FEMData.ndof;    % number of degrees-of-freedom per node
nnp  = FEMData.nnp;     % number of nodal points
nel  = FEMData.nel;     % number of elements
nen  = FEMData.nen;     % number of element nodes

neq  = ndof*nnp;        % Number of equations

f 	= zeros(neq,1);   % Initialize nodal force vector
d 	= zeros(neq,1);   % Initialize nodal displacement matrix
K 	= zeros(neq);     % Initialize stiffness matrix

flags = zeros(neq,1);  % initialize flag vector
e_bc = zeros(neq,1);   % initialize vector of essential boundary condition
n_bc = zeros(neq,1);   % initialize vector of natural boundary condition

% element and material data (given at the element nodes)
E     = FEMData.E;      % nodal values Young's modulus
body  = FEMData.body;   % nodal values body forces
CArea = FEMData.CArea;  % nodal values of cross-sectional area

% gauss integration
ngp = FEMData.ngp;      % number of gauss points

% essential boundary conditions
flags = FEMData.flags;  % flags to mark nodes located on a boundary
                        %   = 2: located on the essential boundary
                        %   = 1: located on the natural boundary
nd = FEMData.nd;        % number of nodes on the essential boundary
e_bc = FEMData.e_bc;    % value of essential B.C
n_bc = FEMData.n_bc;    % value of natural B.C

% point forces
np = FEMData.np;    % number of point forces
if (np>0)
    xp = FEMData.xp;    % array of coordinates where point forces are applied
    P  = FEMData.P;     % array of point forcess
end

% output plots
plot_bar = FEMData.plot_bar;
plot_nod = FEMData.plot_nod;
nplot = nnp*10;    % number of points in the element to plot displacements and stresses

x = FEMData.x;  % x coordinate  
y = FEMData.y;  % y is used only for the bar plot

% connectivity array
IEN = FEMData.IEN;     

plotbar;

% generate LM and ID arrays
d = setup_ID_LM(d);