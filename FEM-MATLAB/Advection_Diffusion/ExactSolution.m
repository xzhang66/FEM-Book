% plots the exact stress  
function ExactSolution 
Global_variables; 

xa = 0:0.01:10;  
ka = k(1);
lene = 0.5;
v = -2*ka*PN/lene;
A = 1/(exp(-10*v/ka)-1);
% plot displacement       
ya = A*exp(-v*xa/ka)-A;
plot(xa,ya,'--r');
