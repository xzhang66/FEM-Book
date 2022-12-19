% plots the exact stress
function ExactSolution 
include_flags; 
 
xx = 0:0.01:2;  
 
subplot(2,1,1); 

% exact displacement for
Ee = 10000;
ue = (-xx.^3/6 + xx)/Ee; 

% plot displacement 
h=plot([xx],[ue], '--r' );       
legend('exact'); 
 
subplot(2,1,2); 

% exact stress
stre = (-xx.^2/2 + 1);

% plot stress 
plot([xx],[stre], '--r' );
