% Gauss-Seidel iteration function. A is matrix of LES, b is the right-side vector
function X = gaussseid(A, b, iterations)

% Determine the number of equations from the length of vector b
N = length(b);

% Initialize the solution vector for the current iteration as zeros
currentIteration = zeros(N, 1);

% Check if the matrix is square
[R, C] = size(A);
if R ~= C
    error("Matrix A must be square.");
end

% Check if the matrix A is diagonally dominant
% This is necessary for ensuring convergence of the method
for m = 1:R
    row = abs(A(m, :));  % Take absolute values of row elements
    diagonal = row(m); % Extract the diagonal element
    offDiagonalSum = sum(row) - diagonal;  % Sum of off-diagonal elements
    if diagonal < offDiagonalSum
        error("Matrix A must be diagonally dominant.");
    end
end

% Perform the Gauss-Seidel iterations
for j = 1:iterations
    for i = 1:N
        % Calculate the sum excluding the diagonal element
        res = sum(A(i, 1:i-1) * currentIteration(1:i-1)) + sum(A(i, i+1:N) * currentIteration(i+1:N));

        % Update the current iteration's solution for i-th element
        currentIteration(i) = (b(i) - res) / A(i, i);
    end

    % Displaying current iteration number
    disp("This is the iteration number: " + j);
end
% Assign the final solution to X
X = currentIteration;
end

% TEST:
% Test with Diagonally Dominant Matrix
%A = [5, -2, 3; -3, 9, 1; 2, -1, 7];
% b = [10; 7; 5];
% iterations = 15;
% X = gaussseid(A, b, iterations);
% disp('Solution vector X:');
% disp(X);

% Test with Non-Diagonally Dominant Matrix
% A = [1, 2, 3; 4, 5, 6; 7, 8, 9];
% b = [14; 32; 50];
% iterations = 10;
% X = gaussseid(A, b, iterations);
% disp('Solution vector X:');
% disp(X);

%  Test with non-square matrix
% A = [1, 2; 3, 4; 5, 6];
% b = [8; 10; 12];
% iterations = 10;
% X = gaussseid(A, b, iterations);
% disp('Solution vector X:');
% disp(X);