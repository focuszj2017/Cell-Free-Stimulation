close all;
clear;

%% Define simulation setup

%Number of setups with random UE locations
nbrOfSetups = 50;

%Number of channel realizations per setup
nbrOfRealizations = 1000;

%Number of APs in the cell-free network
M = 64;

%Number of UEs in the network
K = 16;

%Number of antennas per AP
N = 2;

%Uplink transmit power emmit by each Tag (W)
p = 0.1;

%power control coefficients (W)
alpha_f = 1;

%Carrier frequency(MHz)
f = 900;

%power of noise (W)
segma = 10^(-96/10)/1000;

%Wavelength (m)
lambda = (3*1e8)/(f*1e6);

%Distance between antennas (m)
l = lambda * 0.5;

%Number of APs in the small-cell network
M_cl = 1;

%Number of antennas per AP in small-cell
N_cl = 32;

%Prepare to save simulation results
SE_CF_MRC_tot = zeros(K,nbrOfSetups);
SE_CF_MMSE_tot = zeros(K,nbrOfSetups);
SE_CL_MRC_tot = zeros(K,nbrOfSetups);
SE_CL_MMSE_tot = zeros(K,nbrOfSetups);

%Go through all setups
for n = 1:nbrOfSetups

    %Display simulation progress
    disp(['Setup ' num2str(n) ' out of ' num2str(nbrOfSetups)]);

    %% Cell-Free mMIMO
    %Generate one setup with UEs at random locations
    [Beta,dist,theta] = generateSetup(M,K,N);

    %Generate channel realizations
    %channel gain for MMSE
    [Hhat] = functionComputeChannelGain(nbrOfRealizations,Beta,dist,theta,alpha_f,M,K,N,l,lambda);

    %Compute SE for the Cell-free mMIMO system with Monte Carlo simulations
    [SE_CF_MRC, SE_CF_MMSE] = functionComputeSE(nbrOfRealizations,Beta,dist,theta,Hhat,M,K,N,alpha_f,segma,p,l,lambda);

    %% Small-cell
    %Generate one setup with UEs at random locations
    [Beta,dist,theta] = generateSetup(M_cl,K,N_cl);

    %Generate channel realizations
    %channel gain for MMSE
    [Hhat] = functionComputeChannelGain(nbrOfRealizations, Beta,dist,theta,alpha_f,M_cl,K,N_cl,l,lambda);
    %Compute SE for the small-cell system with Monte Carlo simulations
    [SE_CL_MRC, SE_CL_MMSE] = functionComputeSE(nbrOfRealizations,Beta,dist,theta,Hhat,M_cl,K,N_cl,alpha_f,segma,p,l,lambda);
    
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

% figure(1);
% hold on; box on;
% plot(sort(reshape(SE_CF_MRC_tot,[K*nbrOfSetups,1])), linspace(0,1,K*nbrOfSetups),'r--','LineWidth',2);
% plot(sort(reshape(SE_CF_MMSE_tot,[K*nbrOfSetups,1])), linspace(0,1,K*nbrOfSetups),'r-','LineWidth',2);
% plot(sort(reshape(SE_CL_MRC_tot,[K*nbrOfSetups,1])), linspace(0,1,K*nbrOfSetups),'b--','LineWidth',2);
% plot(sort(reshape(SE_CL_MMSE_tot,[K*nbrOfSetups,1])), linspace(0,1,K*nbrOfSetups),'b-','LineWidth',2);
% xlabel('Spectral Efficiency [bit/s/Hz]');
% ylabel('CDF');
% legend('Cell-Free(MRC)','Cell-Free(MMSE)','Small-Cell(MRC)','Small-Cell(MMSE)');
% 
% figure(2)
% hold on; box on;
% plot(sort(SE_CF_MRC_sum), linspace(0,1,nbrOfSetups),'r--','LineWidth',2);
% plot(sort(SE_CF_MMSE_sum), linspace(0,1,nbrOfSetups),'r-','LineWidth',2);
% plot(sort(SE_CL_MRC_sum), linspace(0,1,nbrOfSetups),'b--','LineWidth',2);
% plot(sort(SE_CL_MMSE_sum), linspace(0,1,nbrOfSetups),'b-','LineWidth',2);
% xlabel('Sum of Spectral Efficency[bit/s/Hz]');
% ylabel('CDF');
% legend('Cell-Free(MRC)','Cell-Free(MMSE)','Small-Cell(MRC)','Small-Cell(MMSE)');

figure(1);
hold on; box on;
plot(sort(reshape(SE_CF_MRC_tot,[K*nbrOfSetups,1])), linspace(0,1,K*nbrOfSetups),'r--','LineWidth',2);
plot(sort(reshape(SE_CF_MMSE_tot,[K*nbrOfSetups,1])), linspace(0,1,K*nbrOfSetups),'r-','LineWidth',2);
plot(sort(reshape(SE_CL_MRC_tot,[K*nbrOfSetups,1])), linspace(0,1,K*nbrOfSetups),'b--','LineWidth',2);
plot(sort(reshape(SE_CL_MMSE_tot,[K*nbrOfSetups,1])), linspace(0,1,K*nbrOfSetups),'b-','LineWidth',2);
xlabel('Spectral Efficiency [bit/s/Hz]');
ylabel('CDF');
legend('Cell-Free(MRC)','Cell-Free(MMSE)','Small-Cell(MRC)','Small-Cell(MMSE)');

figure(2)
hold on; box on;
plot(sort(SE_CF_MRC_sum), linspace(0,1,nbrOfSetups),'r--','LineWidth',2);
plot(sort(SE_CF_MMSE_sum), linspace(0,1,nbrOfSetups),'r-','LineWidth',2);
plot(sort(SE_CL_MRC_sum), linspace(0,1,nbrOfSetups),'b--','LineWidth',2);
plot(sort(SE_CL_MMSE_sum), linspace(0,1,nbrOfSetups),'b-','LineWidth',2);
xlabel('Sum of Spectral Efficency[bit/s/Hz]');
ylabel('CDF');
legend('Cell-Free(MRC)','Cell-Free(MMSE)','Small-Cell(MRC)','Small-Cell(MMSE)');

legend('分布式MIMO接收机（MRC)','分布式MIMO接收机（MMSE）','集中式MIMO接收机（MRC）','集中式MIMO接收机（MMSE）');
ylabel('累计分布函数');
xlabel('频谱效率 [bit/s/Hz]');
xlabel('频谱效率总和 [bit/s/Hz]');
