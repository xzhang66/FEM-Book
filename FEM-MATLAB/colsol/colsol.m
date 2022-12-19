function [IERR, K, R] = colsol(n, m, K, R)
% Solve finite element static equilibrium equations in core using
% the active column solver
%
% Input parameters
%   n      - Number of equations
%   m[n]   - Define the skyline of the stiffness matrix K
%            The row number of the first nonzero element in each column
%   K[n,n] - The stiffness matrix
%   R[n]   - Right-hand-side load vector
%
% Output parameters
%   IERR   - Error indicator. If IERR > 0, K is not positive definite.
%   K      - D and L (Factors of stiffness matrix)
%            The elements of D are stored on its diagonal, 
%            and L replaces its upper triangular part
%   R      - Displacement vector
%
%
% Perform L*D*L(T) factorization of stiffness matrix
%   LDLT is an active column solver to obtain the LDLT factorization
%   of a stiffness matrix K

IERR = 0;
for j = 2:n
    
    for i = m(j)+1:j-1
        c = 0.0;
        for r = max(m(i),m(j)):i-1
            c = c + K(r,i)*K(r,j);
        end
        K(i,j) = K(i,j) - c;
    end
    
    for i = m(j):j-1
        Lij = K(i,j)/K(i,i);
        K(j,j) = K(j,j) - Lij*K(i,j);
        K(i,j) = Lij;
    end

    if K(j,j) <= 0
        fprintf('Error - stiffness matrix is not positive definite !\n')
        fprintf('        Nonpositive pivot for equation %d\n', n)
        fprintf('        Pivot = %f\n', K(j,j))
            
        IERR = j;
        return
    end
end

% Reduce right-hand-side load vector
for i = 2:n
    for j = m(i):i-1
        R(i) = R(i) - K(j,i) * R(j);
    end
end

%back-substitute
for i = 1:n
    R(i) = R(i)/K(i,i);
end

for j = n:-1:2
    for i = m(j):j-1
        R(i) = R(i) - K(i,j)*R(j);
    end
end
