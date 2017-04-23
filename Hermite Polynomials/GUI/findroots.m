function r = findroots(N);
 
J = diag(sqrt([1:N-1]),1)+diag(sqrt([1:N-1]),-1);    % Jacobi matrix
r = eig(sparse(J))/sqrt(2);                    % Compute eigenvalues
