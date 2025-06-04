% Define or load your A, B, and K matrices
% Example: (replace these with your actual system values)
A = [...
     0.00     0.00     0.00     1.00     0.00     0.00     0.00     0.00     0.00;
     0.00     0.00     0.00     0.00     1.00     0.00     0.00     0.00     0.00;
     0.00     0.00     0.00     0.00     0.00     1.00     0.00     0.00     0.00;
   -19.39    -0.95     0.62    -0.47     0.00     0.00     0.00     0.00     0.00;
    -5.87    -5.87     3.87     0.00    -2.89     0.00     0.00     0.00     0.00;
    70.95    70.95    70.95     0.00    -0.00   -53.12     0.00     0.00     0.00;
    -1.00    -0.00    -0.00    -0.00    -0.00    -0.00     0.00     0.00     0.00;
    -0.00    -1.00    -0.00    -0.00    -0.00    -0.00     0.00     0.00     0.00;
    -0.00    -0.00    -1.00    -0.00    -0.00    -0.00     0.00     0.00     0.00];

B = [...
     0.00     0.00     0.00;
     0.00     0.00     0.00;
     0.00     0.00     0.00;
     0.47     0.00     0.00;
     0.00     2.89     0.00;
     0.00     0.00    53.12;
     0.00     0.00     0.00;
     0.00     0.00     0.00;
     0.00     0.00     0.00];

K = [...
   412.01   -1.93    1.32  307.69  -0.01   0.00  -94.28  -0.00  -0.00;
    -1.99  477.77    1.35   -0.02  49.98   0.00    0.00 -47.14  -0.00;
     1.32    1.33  125.71    0.00   0.00  88.04    0.00   0.00 -23.57];

% Compute closed-loop system matrix
A_cl = A - B * K;

% Compute eigenvalues
eigs_cl = eig(A_cl);

% Display results
disp('Closed-loop eigenvalues:');
disp(eigs_cl);

% Optional: plot real vs imaginary
figure;
plot(real(eigs_cl), imag(eigs_cl), 'bx', 'MarkerSize', 10, 'LineWidth', 2);
grid on;
xlabel('Real Part');
ylabel('Imaginary Part');
title('Closed-Loop Poles (Eigenvalues of A - BK)');
