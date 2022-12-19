function plottruss;
include_flags;

% check if truss plot is requested
if strcmpi(plot_truss,'yes')==1;  
    for i = 1:nel
        XX = [x(IEN(i,1)) x(IEN(i,2)) ];
        YY = [y(IEN(i,1)) y(IEN(i,2)) ];
        line(XX,YY);hold on;

        % check if node numbering is requested
        if strcmpi(plot_node,'yes')==1;   
            text(XX(1),YY(1),sprintf('%0.5g',IEN(i,1)));
            text(XX(2),YY(2),sprintf('%0.5g',IEN(i,2)));
        end
    end
    title('Truss Plot');
end

% print mesh parameters
fprintf(1,'\t2D Truss Params \n\n');
fprintf(1,'%s\n\n',Title);
fprintf(1,'No. of Elements  %d \n',nel);
fprintf(1,'No. of Nodes     %d \n',nnp);
fprintf(1,'No. of Equations %d \n\n',neq);