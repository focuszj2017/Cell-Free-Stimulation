close all;
clear;

%% Simulation
%Number of setups with random UE locations
nbrOfSetups = 200;

%Total number of antennas is 36
[SE_CF_MMSE_tot_9_1, SE_CF_MMSE_sum_9_1, SE_CF_MRC_tot_9_1, SE_CF_MRC_sum_9_1] = function_BD_alpha(9,0.1);
[SE_CF_MMSE_tot_9_5, SE_CF_MMSE_sum_9_5, SE_CF_MRC_tot_9_5, SE_CF_MRC_sum_9_5] = function_BD_alpha(9,0.5);
[SE_CF_MMSE_tot_16_1, SE_CF_MMSE_sum_16_1, SE_CF_MRC_tot_16_1, SE_CF_MRC_sum_16_1] = function_BD_alpha(16,0.1);
[SE_CF_MMSE_tot_16_5, SE_CF_MMSE_sum_16_5, SE_CF_MRC_tot_16_5, SE_CF_MRC_sum_16_5] = function_BD_alpha(16,0.5);

%% Plot simulation results

figure(1);
hold on; box on;
plot(sort(reshape(SE_CF_MMSE_tot_9_1,[K*nbrOfSetups,1])), linspace(0,1,9*nbrOfSetups),'r-','LineWidth',2);
plot(sort(reshape(SE_CF_MMSE_tot_9_5,[K*nbrOfSetups,1])), linspace(0,1,9*nbrOfSetups),'b-','LineWidth',2);
plot(sort(reshape(SE_CF_MMSE_tot_16_1,[K*nbrOfSetups,1])), linspace(0,1,16*nbrOfSetups),'m-','LineWidth',2);
plot(sort(reshape(SE_CF_MMSE_tot_16_5,[K*nbrOfSetups,1])), linspace(0,1,16*nbrOfSetups),'k-','LineWidth',2);
xlabel('Spectral Efficiency [bit/s/Hz]');
ylabel('CDF');
legend('K=9,alpha=0.1','K=9,alpha=0.5','K=16,alpha=0.1','K=16,alpha=0.5');

figure(2)
hold on; box on;
plot(sort(SE_CF_MMSE_sum_9_1), linspace(0,1,nbrOfSetups),'r-','LineWidth',2);
plot(sort(SE_CF_MMSE_sum_9_5), linspace(0,1,nbrOfSetups),'b-','LineWidth',2);
plot(sort(SE_CF_MMSE_sum_16_1), linspace(0,1,nbrOfSetups),'m-','LineWidth',2);
plot(sort(SE_CF_MMSE_sum_16_5), linspace(0,1,nbrOfSetups),'k-','LineWidth',2);
xlabel('Sum of Spectral Efficency[bit/s/Hz]');
ylabel('CDF');
legend('K=9,alpha=0.1','K=9,alpha=0.5','K=16,alpha=0.1','K=16,alpha=0.5');