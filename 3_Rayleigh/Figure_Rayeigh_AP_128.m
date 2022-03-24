close all;
clear;

%% Simulation
%Number of setups with random UE locations
nbrOfSetups = 50;

%Number of tags
K = 32;

%Total number of antennas is 128
[SE_CF_MMSE_tot_4_32, SE_CF_MMSE_sum_4_32, SE_CF_MRC_tot_4_32, SE_CF_MRC_sum_4_32] = function_AP_number(4,32,K);
[SE_CF_MMSE_tot_16_8, SE_CF_MMSE_sum_16_8, SE_CF_MRC_tot_16_8, SE_CF_MRC_sum_16_8] = function_AP_number(16,8,K);
[SE_CF_MMSE_tot_64_2, SE_CF_MMSE_sum_64_2, SE_CF_MRC_tot_64_2, SE_CF_MRC_sum_64_2] = function_AP_number(64,2,K);

%% Plot simulation results

figure(1);
hold on; box on;
plot(sort(reshape(SE_CF_MRC_tot_4_32,[K*nbrOfSetups,1])), linspace(0,1,K*nbrOfSetups),'r-','LineWidth',2);
plot(sort(reshape(SE_CF_MRC_tot_16_8,[K*nbrOfSetups,1])), linspace(0,1,K*nbrOfSetups),'b-','LineWidth',2);
plot(sort(reshape(SE_CF_MRC_tot_64_2,[K*nbrOfSetups,1])), linspace(0,1,K*nbrOfSetups),'m-','LineWidth',2);

xlabel('Spectral Efficiency [bit/s/Hz]');
ylabel('CDF');
legend('M=4,N=32','M=16,N=8','M=64,N=2');

figure(2)
hold on; box on;
plot(sort(SE_CF_MRC_sum_4_32), linspace(0,1,nbrOfSetups),'r-','LineWidth',2);
plot(sort(SE_CF_MRC_sum_16_8), linspace(0,1,nbrOfSetups),'b-','LineWidth',2);
plot(sort(SE_CF_MRC_sum_64_2), linspace(0,1,nbrOfSetups),'m-','LineWidth',2);
xlabel('Sum of Spectral Efficency[bit/s/Hz]');
ylabel('CDF');
legend('M=4,N=32','M=16,N=8','M=64,N=2');