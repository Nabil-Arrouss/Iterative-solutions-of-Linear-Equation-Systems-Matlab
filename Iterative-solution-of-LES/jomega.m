% The function examines the parameter of Soothed Jacobi iteration

function [optimalOmega, convergenceInterval, spectralRadii] = jomega(A)
% Initialize range of omega values from 0 to 2 with a step of 0.01
% Omega is the relaxation parameter used in the Soothed Jacobi method.
% convergenceInterval: This is an array or range of omega values for which the Soothed Jacobi method
%                    is expected to converge when applied to the matrix A.
% spectralRadii: This is an array containing the spectral radius calculated for each omega value within
%               the tested range.

omegas = 0:0.01:2;
spectralRadii = zeros(size(omegas));

% Decomposing matrix A into diagonal (D), lower (L), and upper (U) components
 % This decomposition is used in splitting the matrix for the iterative method.
D = diag(diag(A));
L = tril(A, -1);
U = triu(A, 1);

% Iterate over omega values to find the spectral radius
for i = 1:length(omegas)
    omega = omegas(i);

     % Calculate the iteration matrix for the Soothed Jacobi method.
     % This matrix is used to perform the iterative steps in the method.
    M = inv(D + omega * L) * ((1 - omega) * D - omega * U);

    % Compute eigenvalues of the iteration matrix
    eigenvalues = eig(M);

    % Calculate the spectral radius for this omega.
    % The spectral radius is the maximum absolute eigenvalue and is key in determining convergence.
    spectralRadii(i) = max(abs(eigenvalues));
end

% Identify the omega value that minimizes the spectral radius, indicating optimal convergence speed.
[minRadius, optimalIndex] = min(spectralRadii);
optimalOmega = omegas(optimalIndex);

% Determine the interval of omega values that lead to convergence
% Convergence is typically assured when the spectral radius is less than 1.
convergenceInterval = omegas(spectralRadii < 1);

% Plot spectral radius vs omega
figure;
plot(omegas, spectralRadii, 'LineWidth', 2);
hold on;
plot(optimalOmega, minRadius, 'r*', 'MarkerSize', 10);
title('Spectral Radius vs Omega');
xlabel('Omega');
ylabel('Spectral Radius');
grid on;
legend('Spectral Radius', 'Optimal Omega');
hold off;
end


% TEST:
% Test with a Diagonally Dominant Matrix
% A = [4, -1, 0; -1, 4, -1; 0, -1, 3];
% [optimalOmega, convergenceInterval, spectralRadii] = jomega(A);
% disp('Optimal Omega:');
% disp(optimalOmega);
% disp('Convergence Interval:');
% disp(convergenceInterval);

% Test with symmetric positive definite matrix
% A = [5, 2, 1; 2, 6, 3; 1, 3, 4];
% [optimalOmega, convergenceInterval, spectralRadii] = jomega(A);
% disp('Optimal Omega:');
% disp(optimalOmega);
% disp('Convergence Interval:');
% disp(convergenceInterval);

