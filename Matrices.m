A = [0 0 0 1 0 0;
     0 0 0 0 1 0;
     0 0 0 0 0 1;
     -19.3883 -0.9451 0.6223 -0.4659 0 0; 
     -5.8723 -5.8723 3.866 0 -2.8944 0; 
     70.9481 70.9481 70.9481 0 0 -53.1163];

B = [0 0 0;
     0 0 0;
     0 0 0;
     0.4659 0 0;
     0 2.8944 0;
     0 0 53.1163];

K_lqr = [...
  412.01  -1.93   1.32  307.69  -0.01   0.00  -94.28  -0.00  -0.00;
   -1.99 477.77   1.35  -0.02  49.98   0.00   0.00 -47.14  -0.00;
    1.32   1.33 125.71   0.00   0.00  88.04   0.00   0.00 -23.57 ...
];




q_ref = [ 1.75; 1.13; 2.18];