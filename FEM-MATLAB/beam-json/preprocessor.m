% preprocessing: reads input data and sets up mesh information
%   Originally by Haim Waisman, Rensselaer
%   Modified by Xiong Zhang to read input data from json file
%
function  [K,f,d] = preprocessor(DataFile); 
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
neqe = ndof*nen;        % Number of equations for each element 
 
f	= zeros(neq,1);     % Initialize force vector 
d	= zeros(neq,1);     % Initialize displacement vector 
K	= zeros(neq);       % Initialize stiffness matrix 
 
flags = zeros(neq,1);   % initialize flag vector 
e_bc = zeros(neq,1);    % initialize vector of essential boundary condition 
n_bc = zeros(neq,1);    % initialize vector of natural boundary condition 
 
% Element properties  
CArea = FEMData.CArea;  % Elements cross-sectional area   
leng  = FEMData.leng;   % Elements length 
body  = FEMData.body;   % body forces 
E     = FEMData.E;      % Young's Modulus 
 
% gauss integration 
ngp   = FEMData.ngp;    % number of gauss points 
 
% essential boundary conditions 
% odd numbers for displacements; even numbers for rotations 
flags = FEMData.flags;  % flags to mark degrees-of-freedom located on the essential boundary 
e_bc  = FEMData.e_bc;   % value of prescribed displacement/rotation 
nd    = FEMData.nd;     % number of degrees-of-freedom on the essential boundary 
  
% natural boundary conditions 
% odd numbers for shear forces; even numbers for moments 
flags = FEMData.flags;  % flags to mark degrees-of-freedom located on the natural boundary 
n_bc  = FEMData.n_bc;   % value of force 
 
% Applied point forces  
P  = FEMData.P;         % array of point forces           
xp = FEMData.xp;        % array of coordinates where point forces are applied 
np = FEMData.np;        % number of point forces 
 
% output controls 
plot_beam = FEMData.plot_beam; 
plot_nod  = FEMData.plot_nod; 
  
% mesh generation 
x = FEMData.x;     % X coordinate   
y = FEMData.y;     % Y coordinate 
 
% connectivity array 
IEN = FEMData.IEN;      
 
% plot beam 
plotbeam; 

% number of points for plot 
nplot=300; 
 
% Generate LM array  
count = 0; count1 = 0; 
for i = 1:neq 
   if flags(i) == 2                % check if essential boundary 
     count = count + 1; 
     ID(i) = count;                % number first the degrees-of-freedom on essential boundary 
      d(count)= e_bc(i);           % store the reordered values of essential B.C 
   else 
     count1 = count1 + 1; 
     ID(i) = nd + count1; 
   end 
end 
for e = 1:nel 
  for j = 1:nen 
    for m = 1:ndof 
            ind = (j-1)*ndof + m; 
            LM(ind,e) = ID(ndof*IEN(j,e) - ndof + m) ;% create the LM matrix 
    end 
 
  end 
end 
