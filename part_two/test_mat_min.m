%mat = [2 2 2 2 2; 2 1 1 1 2; 2 1 0 1 2; 2 1 1 1 2; 2 2 2 2 2];
%mat = mat+0.1*randn(5);
mat = [7 4 6 5 8
       8 3 2 6 4
       4 1 0 2 5
       5 2 1 4 6
       6 5 8 7 9];
[u_min_index, tau_index] = matrix_minimization(mat);
disp(u_min_index);
disp(tau_index);