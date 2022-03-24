close all;
clear;

%% Simulation
%Number of setups with random UE locations
nbrOfSetups = 2;

%Number of tags
K1 = 9;
K2 = 16;
alpha1 = 0.1;
alpha2 = 0.5;

%Total number of antennas is 36
[SE_CF_MMSE_tot_9_1, SE_CF_MMSE_sum_9_1, SE_CF_MRC_tot_9_1, SE_CF_MRC_sum_9_1] = function_BD_alpha(K1,alpha1);
[SE_CF_MMSE_tot_9_5, SE_CF_MMSE_sum_9_5, SE_CF_MRC_tot_9_5, SE_CF_MRC_sum_9_5] = function_BD_alpha(K1,alpha2);
[SE_CF_MMSE_tot_16_1, SE_CF_MMSE_sum_16_1, SE_CF_MRC_tot_16_1, SE_CF_MRC_sum_16_1] = function_BD_alpha(K2,alpha1);
[SE_CF_MMSE_tot_16_5, SE_CF_MMSE_sum_16_5, SE_CF_MRC_tot_16_5, SE_CF_MRC_sum_16_5] = function_BD_alpha(K2,alpha2);

%% Plot simulation results

figure(1);
hold on; box on;
plot(sort(reshape(SE_CF_MMSE_tot_9_1,[K1*nbrOfSetups,1])), linspace(0,1,K1*nbrOfSetups),'r-','LineWidth',2);
plot(sort(reshape(SE_CF_MMSE_tot_9_5,[K1*nbrOfSetups,1])), linspace(0,1,K1*nbrOfSetups),'b-','LineWidth',2);
plot(sort(reshape(SE_CF_MMSE_tot_16_1,[K2*nbrOfSetups,1])), linspace(0,1,K2*nbrOfSetups),'m-','LineWidth',2);
plot(sort(reshape(SE_CF_MMSE_tot_16_5,[K2*nbrOfSetups,1])), linspace(0,1,K2*nbrOfSetups),'k-','LineWidth',2);
xlabel('Spectral Efficiency [bit/s/Hz]');
ylabel('CDF');
legend(['K=',num2str(K1),',alpha=0.1'],['K=',num2str(K1),',alpha=0.5'],['K=',num2str(K2),',alpha=0.1'],['K=',num2str(K2),',alpha=0.5']);

figure(2)
hold on; box on;
plot(sort(SE_CF_MMSE_sum_9_1), linspace(0,1,nbrOfSetups),'r-','LineWidth',2);
plot(sort(SE_CF_MMSE_sum_9_5), linspace(0,1,nbrOfSetups),'b-','LineWidth',2);
plot(sort(SE_CF_MMSE_sum_16_1), linspace(0,1,nbrOfSetups),'m-','LineWidth',2);
plot(sort(SE_CF_MMSE_sum_16_5), linspace(0,1,nbrOfSetups),'k-','LineWidth',2);
xlabel('Sum of Spectral Efficency[bit/s/Hz]');
ylabel('CDF');
legend(['K= ',num2str(K1),',alpha=0.1'],['K=',num2str(K1),',alpha=0.5'],['K=',num2str(K2),',alpha=0.1'],['K=',num2str(K2),',alpha=0.5']);