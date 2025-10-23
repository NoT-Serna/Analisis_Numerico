function [xx, uniqueSol] = iterativeRefinement(A, b, x0, N)
    % iterativeRefinement performs N steps of iterative refinement
    % starting from an initial solution x0.
    %
    % Inputs:
    %   A  - Coefficient matrix (n x n)
    %   b  - Right-hand side vector (n x 1)
    %   x0 - Initial solution (n x 1)
    %   N  - Number of refinement iterations
    %
    % Outputs:
    %   xx        - Refined solution after N iterations
    %   uniqueSol - Logical flag indicating if the solution converged

    n = size(A, 1);
    xx = x0;          % Start from given initial guess
    uniqueSol = true;

    % LU factorization (done once)
    [L, U, P] = lu(A);

    % Iterative refinement loop
    for k = 1:N
        r = b - A * xx;            % Residual
        y = U \ (L \ (P * r));     % Correction using LU
        xx = xx + y;               % Update solution
    end

    % Check if the residual is small enough
    if norm(b - A * xx, inf) > 1e-8
        uniqueSol = false;
        disp('Iterative refinement finished, but residual remains high.');
    end
end
