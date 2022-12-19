%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2D Elasticity matlab code from chapter 9                       %
%    J. Fish & T. Belytschko: A First Course in Finite Elements  %
% Programed by: Haim Waisman, Rensselaer                         %
% Revised by:   Xiong Zhang to read data from json file          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; 
close all;  
 
% include global variables 
include_flags; 
 
% Preprocessing  K
[K,f,d] = preprocessor('elasticity_16.json'); 
 
% Calculation and assembly of element matrices
for e = 1:nel 
    [ke, fe] = elast2Delem(e);  
    [K,f] = assembly(K,f,e,ke,fe); 
end 
 
% Compute and assemble nodal boundary force vector and point forces 
f = point_and_trac(f); 
 
% Solution Phase 
[d,f_E] = solvedr(K,f,d); 
 
% Postprocessor 
postprocess(d);
