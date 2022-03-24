close all;
clear;

%% Simulation
%Number of setups with random UE locations
nbrOfSetups = 100;

%Number of tags
K = 9;

%Total number of antennas is 36
[SE_CF_MMSE_tot_4_9, SE_CF_MMSE_sum_4_9, SE_CF_MRC_tot_4_9, SE_CF_MRC_sum_4_9] = function_AP_number(4,9,K,nbrOfSetups);
[SE_CF_MMSE_tot_9_4, SE_CF_MMSE_sum_9_4, SE_CF_MRC_tot_9_4, SE_CF_MRC_sum_9_4] = function_AP_number(9,4,K,nbrOfSetups);
[SE_CF_MMSE_tot_1_36, SE_CF_MMSE_sum_1_36, SE_CF_MRC_tot_1_36, SE_CF_MRC_sum_1_36] = function_AP_number(1,36,K,nbrOfSetups);
[SE_CF_MMSE_tot_36_1, SE_CF_MMSE_sum_36_1, SE_CF_MRC_tot_36_1, SE_CF_MRC_sum_36_1] = function_AP_number(36,1,K,nbrOfSetups);

%% Plot simulation results

figure(1);
hold on; box on;
plot(sort(reshape(SE_CF_MMSE_tot_1_36,[K*nbrOfSetups,1])), linspace(0,1,K*nbrOfSetups),'r-','LineWidth',2);
plot(sort(reshape(SE_CF_MMSE_tot_4_9,[K*nbrOfSetups,1])), linspace(0,1,K*nbrOfSetups),'b-','LineWidth',2);
plot(sort(reshape(SE_CF_MMSE_tot_9_4,[K*nbrOfSetups,1])), linspace(0,1,K*nbrOfSetups),'m-','LineWidth',2);
plot(sort(reshape(SE_CF_MMSE_tot_36_1,[K*nbrOfSetups,1])), linspace(0,1,K*nbrOfSetups),'k-','LineWidth',2);
xlabel('Spectral Efficiency [bit/s/Hz]');
ylabel('CDF');
legend('M=1,N=36','M=4,N=9','M=9,N=4','M=36,N=1');

figure(2)
hold on; box on;
plot(sort(SE_CF_MMSE_sum_1_36), linspace(0,1,nbrOfSetups),'r-','LineWidth',2);
plot(sort(SE_CF_MMSE_sum_4_9), linspace(0,1,nbrOfSetups),'b-','LineWidth',2);
plot(sort(SE_CF_MMSE_sum_9_4), linspace(0,1,nbrOfSetups),'m-','LineWidth',2);
plot(sort(SE_CF_MMSE_sum_36_1), linspace(0,1,nbrOfSetups),'k-','LineWidth',2);
xlabel('Sum of Spectral Efficency[bit/s/Hz]');
ylabel('CDF');
legend('M=1,N=36','M=4,N=9','M=9,N=4','M=36,N=1');