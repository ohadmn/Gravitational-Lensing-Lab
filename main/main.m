% main.m
% Define the file path to the data file
filePath = "~/Desktop/file_one.dat";
 

% Call the processing function and retrieve results
% [max_time, flux_max] = leastSquaresPolynomialFit(filePath);

% part 2 section d.(Calculate the mass from eq.1.13 v=150 km/s ; x = 1/2 ; Ds = 8kpc)
% disp(filePath)
% [tau_final_file_one] = deltaTau(filePath, 2454476.072, 51.032, 0.121);
% disp(tau_final_file_one);
% 
% filePath = "~/Desktop/BLG012.dat";
% disp(filePath)
% [tau_final_BGL012] = deltaTau(filePath, 2454501.089, 9.322, 0.576);
% disp(tau_final_BGL012);
% 
% filePath = "~/Desktop/BLG004.dat";
% disp(filePath)
% [tau_final_BGL004] = deltaTau(filePath, 2454521.168, 18.108, 0.145);
% disp(tau_final_BGL004);
% 
% filePath = "~/Desktop/BLG005.dat";
% disp(filePath)
% [tau_final_BGL005] = deltaTau(filePath, 2454534.171, 17.577, 0.505);
% disp(tau_final_BGL005);
% 
% filePath = "~/Desktop/BLG007.dat";
% disp(filePath)
% [tau_final_BGL007] = deltaTau(filePath, 2454534.013, 28.01, 1.397);
% disp(tau_final_BGL007);
% 
% filePath = "~/Desktop/BLG009.dat";
% disp(filePath)
% [tau_final_BLG009] = deltaTau(filePath, 2454529.749, 28.962, 0.788);
% disp(tau_final_BLG009);


% part 3
% [tau_final_file_one, u_min_final_file_one, T0_final_file_one] = part3(filePath, 2454476.072, 51.032, 0.121);
% disp(tau_final_file_one);
% disp(u_min_final_file_one);
% disp(T0_final_file_one);

%part 3 test
% [tau_final_file_one, u_min_final_file_one, T0_final_file_one] = Copy_of_part3(filePath, 2454476.072, 51.032, 0.121);
% disp(tau_final_file_one);
% disp(u_min_final_file_one);
% disp(T0_final_file_one);

% part 4
% filePath = "~/Desktop/BLG002.dat";
% [tau_final_file_one, u_min_final_file_one, fbl_final_file_one] = part4(filePath, 2454517.540, 0.512, 16.779, 0.411);
% disp(tau_final_file_one);
% disp(u_min_final_file_one);
% disp(fbl_final_file_one);

% test part 4
% filePath = "~/Desktop/BLG002.dat";
% [tau_final_file_one, u_min_final_file_one, fbl_final_file_one] = part4(filePath, 2454517.540, 0.512, 16.779, 0.411);
% disp(tau_final_file_one);
% disp(u_min_final_file_one);
% disp(fbl_final_file_one);

% part 5
% filePath = "~/Desktop/BLG002.dat";
% [tau_final_file_one, u_min_final_file_one, fbl_final_file_one] = part5(filePath, 2454517.540, 0.512, 16.779, 0.411);
% disp(tau_final_file_one);
% disp(u_min_final_file_one);
% disp(fbl_final_file_one);

% test part 5
filePath = "~/Desktop/BLG002.dat";
[tau_final_file_one, u_min_final_file_one, fbl_final_file_one, T0_final_file_one] = Copy_of_part5(filePath, 2454517.540, 0.512, 16.779, 0.411);
disp(tau_final_file_one);
disp(u_min_final_file_one);
disp(fbl_final_file_one);
disp(T0_final_file_one);