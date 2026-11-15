% % main.m
% % Define the file path to the data file
filePath = "~/Desktop/SMC001.dat";
% data = load(filePath);
% 
% % part 1
% % Call the processing function and retrieve results
% [max_time, flux_max] = leastSquaresPolynomialFit(filePath);
% 
% % part 2
% [tau_final_SMC001, u_min_final_SMC001] = contourErrorGraph(filePath, 2454476.072, 51.032, 0.121);
% disp(tau_final_SMC001);
% disp(['tau_SMC001 accuracy:', num2str((tau_final_SMC001/51.032)*100)]);
% disp(u_min_final_SMC001);
% disp(['u_min_SMC001 accuracy:', num2str((u_min_final_SMC001/0.121)*100)]);
% 
% % part 2 section d.(Calculate the mass from eq.1.13 v=150 km/s ; x = 1/2 ; Ds = 8kpc)
% disp(filePath)
% [tau_final_SMC001] = deltaTau(filePath, 2454476.072, 51.032, 0.121);
% disp(['tau_SMC001 accuracy:', num2str((tau_final_SMC001/51.032)*100)]);
% 
% filePath = "~/Desktop/BLG004.dat";
% disp(filePath)
% [tau_final_BLG004] = deltaTau(filePath, 2454521.168, 18.108, 0.145);
% disp(['tau_BLG004 accuracy:', num2str((tau_final_BLG004/18.108)*100)]);
% 
% filePath = "~/Desktop/BLG005.dat";
% disp(filePath)
% [tau_final_BLG005] = deltaTau(filePath, 2454534.171, 17.577, 0.505);
% disp(['tau_BLG005 accuracy:', num2str((tau_final_BLG005/17.577)*100)]);
% 
% filePath = "~/Desktop/BLG007.dat";
% disp(filePath)
% [tau_final_BLG007] = deltaTau(filePath, 2454534.013, 28.01, 1.397);
% disp(['tau_BLG007 accuracy:', num2str((tau_final_BLG007/28.01)*100)]);
% 
% filePath = "~/Desktop/BLG012.dat";
% disp(filePath)
% [tau_final_BLG012] = deltaTau(filePath, 2454501.089, 9.322, 0.576);
% disp(['tau_BLG012 accuracy:', num2str((tau_final_BLG012/9.322)*100)]);
% 
% filePath = "~/Desktop/BLG015.dat";
% disp(filePath)
% [tau_final_BLG015] = deltaTau(filePath, 2454510.512, 6.520, 1.105);
% disp(['tau_BLG015 accuracy:', num2str((tau_final_BLG015/6.520)*100)]);
% 
% 
% % part 3 
% filePath = "~/Desktop/SMC001.dat";
% [tau_final_SMC001, u_min_final_SMC001, T0_final_SMC001] = part3Vectorized(filePath, 2454476.072, 51.032, 0.121);
% disp(tau_final_SMC001);
% disp(['tau_SMC001 accuracy:', num2str((tau_final_SMC001/51.032)*100)]);
% disp(u_min_final_SMC001);
% disp(['u_min_SMC001 accuracy:', num2str((u_min_final_SMC001/0.121)*100)]);
% disp(T0_final_SMC001);
% disp(['T0_SMC001 accuracy:', num2str((T0_final_SMC001/2454476.072)*100)]);
% 
% % test part 4
% filePath = "~/Desktop/BLG002.dat";
% [tau_final_BLG002, u_min_final_BLG002, fbl_final_BLG002]  = vector_part4(filePath, 2454517.540, 0.512, 16.779, 0.411);
% disp(tau_final_BLG002);
% disp(['tau_BLG0002 accuracy:', num2str((tau_final_BLG002/16.779)*100)]);
% disp(u_min_final_BLG002);
% disp(['u_min_BLG0002 accuracy:', num2str((u_min_final_BLG002/0.411)*100)]);
% disp(fbl_final_BLG002);
% disp(['fbl_BLG0002 accuracy:', num2str((fbl_final_BLG002/0.512)*100)]);
% 
% % test part 5
% filePath = "~/Desktop/BLG002.dat";
% [tau_final_BLG002, u_min_final_BLG002, fbl_final_BLG002, T0_final_BLG002] = vector_part5(filePath, 2454517.540, 0.512, 16.779, 0.411);
% disp(tau_final_BLG002);
% disp(['tau_BLG0002 accuracy:', num2str((tau_final_BLG002/16.779)*100)]);
% disp(u_min_final_BLG002);
% disp(['u_min_BLG0002 accuracy:', num2str((u_min_final_BLG002/0.411)*100)]);
% disp(fbl_final_BLG002);
% disp(['fbl_BLG0002 accuracy:', num2str((fbl_final_BLG002/0.512)*100)]);
% disp(T0_final_BLG002);
% disp(['T0_BLG0002 accuracy:', num2str((T0_final_BLG002/2454476.072)*100)]);
% 
% filePath = "~/Desktop/BLG002.dat";
% part5_vec_contour_slice(filePath, T0_final_BLG002, fbl_final_BLG002, tau_final_BLG002, u_min_final_BLG002);
% 
% test part 6
filePath = "~/Desktop/BLG002.dat";
[tau_final_BLG002, u_min_final_BLG002, fbl_final_BLG002, T0_final_BLG002] = part6(filePath, 2454517.540, 0.512, 16.779, 0.411);
disp(tau_final_BLG002);
disp(['tau_BLG0002 accuracy:', num2str((tau_final_BLG002/16.779)*100)]);
disp(u_min_final_BLG002);
disp(['u_min_BLG0002 accuracy:', num2str((u_min_final_BLG002/0.411)*100)]);
disp(fbl_final_BLG002);
disp(['fbl_BLG0002 accuracy:', num2str((fbl_final_BLG002/0.512)*100)]);
disp(T0_final_BLG002);
disp(['T0_BLG0002 accuracy:', num2str((T0_final_BLG002/2454476.072)*100)]);
% 
% 
% 
% % NOT_WORKING
% % filePath = "~/Desktop/BLG009.dat";
% % disp(filePath)
% % [tau_final_BLG009] = deltaTau(filePath, 2454529.749, 28.962, 0.788);
% % filePath = "~/Desktop/BLG016.dat";
% % disp(filePath)
% % [tau_final_BLG016] = deltaTau(filePath, 2454515.183, 6.371, 0.332);
% 
% % % Archive 
% % % part 5
% % filePath = "~/Desktop/BLG002.dat";
% % [tau_final_file_one, u_min_final_file_one, fbl_final_file_one] = part5(filePath, 2454517.540, 0.512, 16.779, 0.411);
% % disp(tau_final_file_one);
% % disp(u_min_final_file_one);
% % disp(fbl_final_file_one);
% 
% 
% % part 3
% % [tau_final_file_one, u_min_final_file_one, T0_final_file_one] = part3(filePath, 2454476.072, 51.032, 0.121);
% % disp(tau_final_file_one);
% % disp(u_min_final_file_one);
% % disp(T0_final_file_one);
% % 
% % part 4
% % filePath = "~/Desktop/BLG002.dat";
% % [tau_final_BLG002, u_min_final_BLG002, fbl_final_BLG002] = part4(filePath, 2454517.540, 0.512, 16.779, 0.411);
% % disp(tau_final_BLG002);
% % disp(u_min_final_BLG002);
% % disp(fbl_final_BLG002);
% % 
% % filePath = "~/Desktop/BLG002.dat";
% % [tau_final_file_one, u_min_final_file_one, fbl_final_file_one, T0_final_file_one] = Copy_of_part6(filePath, 2454517.540, 0.512, 16.779, 0.411);
% % disp(tau_final_file_one);
% % disp(u_min_final_file_one);
% % disp(fbl_final_file_one);
% % disp(T0_final_file_one);
% 
% % filePath = "~/Desktop/BLG016.dat";
% % disp(filePath)
% % [tau_final_BLG016] = deltaTau(filePath, 2454515.183, 6.371, 0.332);
% % disp(['tau_BLG016 accuracy:', num2str((tau_final_BLG016/6.371)*100)]);