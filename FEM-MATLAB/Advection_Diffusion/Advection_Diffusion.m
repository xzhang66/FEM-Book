%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                %
% 1D FEM Program for Advection-Diffusion Problem %
% with Galerkin and Petrov-Galerkin Method       %
% ( Modified from Program bar1D)                 %
%                                                %
% Computational Dynamics Lab                     %
% 2019.3.19                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear all;

% Include global variables
Global_variables;

% Preprocessing
[K,f,d] = preprocessor;

dd = zeros(neq,3);
xx = zeros(neq,1);

for i=1:size(alpha,2)

% Element matrix computations and assembly
K=zeros(neq);
for e = 1:nel
    [ke,fe] = element(e,alpha(i));
    [K, f] = assembly(K,f,e,ke,fe);
end

% Add nodal boundary force vector
f = naturalBC(f);

% Partition and solution
[d,f_E] = solvedr(K,f,d);

% Postprocessing
postprocessor(d);

for j=1:nel
    dd(j,i) = d(IEN(1,j));
end
dd(neq,i) = d(IEN(2,nel));

end

for j=1:nel
    xx(j) = x(IEN(1,j));
end
xx(neq) = x(IEN(2,nel));

fname = ['AD-Peclet-' num2str(PN) '.dat'];
fid = fopen(fname,'w');

for i=1:neq
    fprintf(fid,'%f  %f  %f  %f\n', xx(i),dd(i,1), dd(i,2), dd(i,3));
end

fclose(fid);