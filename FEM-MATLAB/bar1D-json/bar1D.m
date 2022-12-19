%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1D FEM Program (Chapter 5)                                   %
%   Originally by: Haim Waisman, Rensselaer                    %
%   Modified by: Xiong Zhang to read input data from json file %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;

% include global variables
include_flags;

% Preprocessing
[K,f,d] = preprocessor('bar_5_2_2.json');

% Preporcessing for convergence analysis
% [K,f,d] = preprocessor('convergence/16-elements-3Q.json');

% Element matrix computations and assembly
for e = 1:nel
    [ke,fe] = barelem(e);
    [K, f] = assembly(K,f,e,ke,fe);
end

% Add nodal boundary force vector
f = naturalBC(f);

% Partition and solution
[d,f_E] = solvedr(K,f,d);

% Postprocessing
postprocessor(d);

% plot the exact solution
ExactSolution;

% Convergence study
% ExactSolution_convergence;
% ErrorNorm(d);
