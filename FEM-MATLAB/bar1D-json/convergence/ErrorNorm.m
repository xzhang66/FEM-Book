% Calculate Error norm (L2 and energy norm)
function ErrorNorm(d)
include_flags;

ngp = 3;
[w, gp] = gauss(ngp);    % extract Gauss points and weights

L2Norm = 0;
EnNorm = 0;

L2NormEx = 0;
EnNormEx = 0;

for e = 1:nel
    
    de = d(LM(:,e));          % extract element nodal displacements
    IENe = IEN(:,e);          % extract local connectivity information
    xe = x(IENe);             % extract element x coordinates
    J = (xe(nen) - xe(1))/2;  % compute Jacobian
    
    for i = 1:ngp
        xt = 0.5*(xe(1)+xe(nen))+J*gp(i);  % Compute Gauss points in physical coordinates
        
        N = Nmatrix1D(xt,xe);    % shape functions matrix
        B = Bmatrix1D(xt,xe);    % derivative of shape functions matrix
        
        Ee = N*E(IENe);          % Young's modulus at element gauss points
        
        uh  = N*de;              % displacement at gauss point
        uex = (-xt^3/6 + xt)/Ee; % Exact displacement
        L2Norm = L2Norm + J*w(i)*(uex - uh)^2;
        L2NormEx = L2NormEx + J*w(i)*(uex)^2;
        
        sh  = B*de;              % strain at Gauss points
        sex = (-xt^2/2 + 1)/Ee;  % Exact strain
        EnNorm = EnNorm + 0.5*J*w(i)*Ee*(sex-sh)^2;
        EnNormEx = EnNormEx + 0.5*J*w(i)*Ee*(sex)^2;
    end
    
end

L2Norm = sqrt(L2Norm);
L2NormEx = sqrt(L2NormEx);

EnNorm = sqrt(EnNorm);
EnNormEx = sqrt(EnNormEx);

% print stresses at element gauss points
fprintf(1,'\nErro norms\n');
fprintf(1,'h\tL2Norm\t\tL2NormRel\tEnNorm\t\tEnNormRel\n');
fprintf(1,'%g\t%g\t%g\t%g\t%g\n',2/nel,L2Norm,L2Norm/L2NormEx,EnNorm,EnNorm/EnNormEx);
