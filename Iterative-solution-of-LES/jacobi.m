% Jacobi iteration function. A is matrix of LES, b is the right-side vector
function X = jacobi(A,b,iterations)

% X is the approx. of the solution vector

% To determine number of equations by the lenghth of vector b
N = length(b);

% Initialize the solution vector for the previous iteration as zeros
previousIteration = zeros(N,1);

% Initialize the solution vector for the current iteration as zeros
currentIteration = zeros(N,1);

% CHeck if the matrix is square!
[R,C] = size(A);

if R ~= C
    error("Matrix A must be square.");
end

% Check if the matrix A is diagonally dominant
for m = 1:R
    row = abs(A(m, :)); % Take absolute values of row elements
    diagonal = row(m);  % Extract the diagonal element
    offDiagonalSum = sum(row) - diagonal;  % Sum of off-diagonal elements
    if diagonal < offDiagonalSum
        error('Matrix A must be diagonally dominant.');
    end
end

% Performing the Jacobi iterations
for j = 1:iterations
    for i = 1:N
        % Calculate the sum excluding the diagonal element
        res = sum(A(i, [1:i-1, i+1:N]) * previousIteration([1:i-1, i+1:N]));

        % Update the current iteration's solution for i-th element
        currentIteration(i) = (b(i) - res) / A(i, i);
    end

    % Displaying current iteration number
    disp("The current iteration number is: " + j);

    % Updating solution vector for the next iteration
    previousIteration=currentIteration;
end
 % Return the final iteration result as the solution
X = currentIteration;
end

% TEST:
% Test with diagonally dominant matrix
% A = [4, -1, 0, 0; -1, 4, -1, 0; 0, -1, 4, -1; 0, 0, -1, 3];
% b = [15; 10; 10; 10];
% iterations = 25;
% X = jacobi(A, b, iterations);
% disp('Solution vector X:');
% disp(X);

% Test with non-diagonally dominant matrix
%  A = [1, 2, 3; 4, 5, 6; 7, 8, 9];
% b = [14; 32; 50];
% iterations = 10;
% X = jacobi(A, b, iterations);
% disp('Solution vector X:');
% disp(X);

% Test non-square matrix
%  A = [1, 2; 3, 4; 5, 6];
% b = [8; 10; 12];
% iterations = 10;
% X = jacobi(A, b, iterations);
% disp('Solution vector X:');
% disp(X);