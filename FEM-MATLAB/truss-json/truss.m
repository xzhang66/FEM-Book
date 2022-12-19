%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2D Truss matlab code from                                      %
%    J. Fish & T. Belytschko: A First Course in Finite Elements  %
% Programed by: Haim Waisman, Rensselaer                         %
% Revised by:   Xiong Zhang to read data from json file          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all; 
 
% include global variables
include_flags;  

% Preprocessor Phase 
 [K,f,d]	= preprocessor('truss_2_1.json');

% Calculation and assembly of element matrices
for e = 1:nel
    ke	= trusselem(e);
    K	= assembly(K,e,ke);
end

% Solution Phase
 [d,f_E]	= solvedr(K,f,d);
 
% Postprocessor Phase 
postprocessor(d)