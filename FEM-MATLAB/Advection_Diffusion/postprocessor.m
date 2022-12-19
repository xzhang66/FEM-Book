% postprocessing
function postprocessor(d)
Global_variables;


figure();
hold on;

for i=1:nel
    plot(x(IEN(:,i)),d(IEN(:,i)),'-*');
end
ExactSolution;

end
