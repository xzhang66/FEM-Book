% Read in equations from JSON file
Equations = jsondecode(fileread('Example_n_5.json'));

n = Equations.n;    % Number of equations
m = Equations.m;    % The row number of the first nonzero element in each column
K = Equations.K;    % The stiffness matrix
R = Equations.R;    % The Right-hand-side load vector

[IERR, K, R] = colsol(n, m, K, R);