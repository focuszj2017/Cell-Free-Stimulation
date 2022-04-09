close all;
clear;

%% Define simulation setup

%Number of setups with random UE locations
nbrOfSetups = 200;

%Number of channel realizations per setup
nbrOfRealizations = 2000;

%Number of APs in the cell-free network
M = 16;

%Number of UEs in the network
K = 4;

%Number of antennas per AP
N = 2;

%Uplink transmit power emmit by each Tag (W)
p = 0.1;

%power control coefficients (W)
alpha_f = 1;

%power of noise (W)
segma = 10^(-96/10)/1000;

%Number of APs in the small-cell network
M_cl = 1;

%Number of antennas per AP in small-cell
N_cl = 32;

%Prepare to save simulation results
SE_CF_MRC_tot = zeros(K,nbrOfSetups);
SE_CF_MMSE_tot = zeros(K,nbrOfSetups);
SE_CL_MRC_tot = zeros(K,nbrOfSetups);
SE_CL_MMSE_tot = zeros(K,nbrOfSetups);

SE_CF_MRC_sum_tot = zeros(nbrOfRealizations,nbrOfSetups);
SE_CF_MMSE_sum_tot = zeros(nbrOfRealizations,nbrOfSetups);
SE_CL_MRC_sum_tot = zeros(nbrOfRealizations,nbrOfSetups);
SE_CL_MMSE_sum_tot = zeros(nbrOfRealizations,nbrOfSetups);

%Go through all setups
for n = 1:nbrOfSetups

    %Display simulation progress
    %disp(['Setup ' num2str(n) ' out of ' num2str(nbrOfSetups)]);

    %% Cell-Free mMIMO
    %Generate one setup with UEs at random locations
    [pathLoss] = generateSetup(M,K,N);

    %Generate channel realizations
    %channel gain for MMSE
    [Hhat] = functionComputeChannelGain(nbrOfRealizations, pathLoss, alpha_f, M, K, N);

    %Compute SE for the Cell-free mMIMO system with Monte Carlo simulations
    [SE_CF_MRC, SE_CF_MMSE] = functionComputeSE_parfor(nbrOfRealizations,pathLoss,Hhat,M,K,N,alpha_f,segma,p);

    %% Small-cell
    %Generate one setup with UEs at random locations
    [pathLoss] = generateSetup(M_cl,K,N_cl);

    %Generate channel realizations
    %channel gain for MMSE
    [Hhat] = functionComputeChannelGain(nbrOfRealizations, pathLoss, alpha_f, M_cl, K, N_cl);
    %Compute SE for the small-cell system with Monte Carlo simulations
    [SE_CL_MRC, SE_CL_MMSE] = functionComputeSE_parfor(nbrOfRealizations,pathLoss,Hhat,M_cl,K,N_cl,alpha_f,segma,p);
    
    %Save SE values
    SE_CF_MRC_tot(:,n) = SE_CF_MRC;
    SE_CF_MMSE_tot(:,n) = SE_CF_MMSE;
    SE_CL_MRC_tot(:,n) = SE_CL_MRC;
    SE_CL_MMSE_tot(:,n) = SE_CL_MMSE;
    
    %Remove large matrices at the end of analyzing this setup
    %clear Hhat;

end

SE_CF_MRC_sum = sum(SE_CF_MRC_tot,1);
SE_CF_MMSE_sum = sum(SE_CF_MMSE_tot,1);
SE_CL_MRC_sum = sum(SE_CL_MRC_tot,1);
SE_CL_MMSE_sum = sum(SE_CL_MMSE_tot,1);

%% Plot simulation results

figure(1);
hold on; box on;
plot(sort(reshape(SE_CF_MRC_tot,[K*nbrOfSetups,1])), linspace(0,1,K*nbrOfSetups),'r--','LineWidth',2);
plot(sort(reshape(SE_CF_MMSE_tot,[K*nbrOfSetups,1])), linspace(0,1,K*nbrOfSetups),'r-','LineWidth',2);
plot(sort(reshape(SE_CL_MRC_tot,[K*nbrOfSetups,1])), linspace(0,1,K*nbrOfSetups),'b--','LineWidth',2);
plot(sort(reshape(SE_CL_MMSE_tot,[K*nbrOfSetups,1])), linspace(0,1,K*nbrOfSetups),'b-','LineWidth',2);
xlabel('频谱效率[bit/s/Hz]');
ylabel('累计分布函数');
legend('分布式MIMO(MRC)','分布式MIMO(MMSE)','集中式MIMO(MRC)','集中式MIMO(MMSE)');

figure(2)
hold on; box on;
plot(sort(SE_CF_MRC_sum), linspace(0,1,nbrOfSetups),'r--','LineWidth',2);
plot(sort(SE_CF_MMSE_sum), linspace(0,1,nbrOfSetups),'r-','LineWidth',2);
plot(sort(SE_CL_MRC_sum), linspace(0,1,nbrOfSetups),'b--','LineWidth',2);
plot(sort(SE_CL_MMSE_sum), linspace(0,1,nbrOfSetups),'b-','LineWidth',2);
xlabel('频谱效率总和[bit/s/Hz]');
ylabel('累计分布函数');
legend('分布式MIMO(MRC)','分布式MIMO(MMSE)','集中式MIMO(MRC)','集中式MIMO(MMSE)');


legend('分布式MIMO接收机（MRC)','分布式MIMO接收机（MMSE）','集中式MIMO接收机（MRC）','集中式MIMO接收机（MMSE）');
ylabel('累计分布函数');
xlabel('频谱效率 [bit/s/Hz]');
xlabel('频谱效率总和 [bit/s/Hz]');

legend('分布式MIMO系统','集中式MIMO系统')